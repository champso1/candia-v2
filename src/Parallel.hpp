/**
 *  @file Parallel.hpp
 *  @brief Internal helpers for running work in parallel while propagating exceptions to the caller.
 *
 *  An exception that escapes a std::thread or a std::execution::par algorithm calls std::terminate,
 *  so these helpers catch the first exception thrown by any worker and rethrow it once all work is done.
 */

#pragma once

#include <algorithm>
#include <exception>
#include <execution>
#include <functional>
#include <mutex>
#include <ranges>
#include <thread>
#include <utility>
#include <vector>

namespace Candia2
{
	/** @brief Stores the first exception reported by any worker. */
	class FirstException final
	{
		std::exception_ptr _error{};
		std::mutex _mutex{};
	public:
		/** @brief records the exception currently being handled, unless one was already recorded */
		void capture()
		{
			std::lock_guard<std::mutex> guard{_mutex};
			if (!_error)
				_error = std::current_exception();
		}
		/** @brief rethrows the recorded exception, if any */
		void rethrow()
		{
			if (_error)
				std::rethrow_exception(std::exchange(_error, nullptr));
		}
	};

	/**
	 *  @brief Calls @a body for every element of @a range using std::execution::par
	 *  @note rethrows the first exception thrown by @a body after every call has finished
	 */
	template <std::ranges::range TRange, typename TBody>
	void parallelFor(TRange&& range, TBody&& body)
	{
		FirstException error{};
		std::for_each(std::execution::par, std::ranges::begin(range), std::ranges::end(range), [&](auto&& elem) {
			try {
				body(elem);
			} catch (...) {
				error.capture();
			}
		});
		error.rethrow();
	}

	/**
	 *  @brief A group of std::threads that rethrows the first exception from any of them on @a join()
	 */
	class ThreadGroup final
	{
		std::vector<std::thread> _threads{};
		FirstException _error{};
	public:
		ThreadGroup() = default;
		ThreadGroup(ThreadGroup const&) = delete;
		ThreadGroup& operator=(ThreadGroup const&) = delete;
		~ThreadGroup()
		{
			for (std::thread& t : _threads)
				if (t.joinable())
					t.join();
		}

		/** @brief starts a thread running std::invoke(func, args...) */
		template <typename TFunc, typename... TArgs>
		void spawn(TFunc func, TArgs... args)
		{
			_threads.emplace_back([this, func, args...]() {
				try {
					std::invoke(func, args...);
				} catch (...) {
					_error.capture();
				}
			});
		}

		/** @brief waits for every thread, then rethrows the first exception any of them threw */
		void join()
		{
			for (std::thread& t : _threads)
				t.join();
			_threads.clear();
			_error.rethrow();
		}
	};
}
