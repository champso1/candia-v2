/**
 *  @file Common.hpp
 *  @brief Contains type definitions, constants, logging, and other misc definitions.
 */

#ifndef __COMMON_HPP
#define __COMMON_HPP

#include <functional>
#include <vector>
#include <numbers>
#include <format>
#include <iostream>
#include <sstream>
#include <iterator>

using uint = unsigned;

#define UNUSED(x) (void)(x)
#define EPS 1e-8


namespace Candia2
{
	/**
	 *  @defgroup constants Constants
	 *  @{
	 */
	constexpr const double CF    = 4.0/3.0;
	constexpr const double NC    = 3.0;
	constexpr const double TR    = 0.5;
	constexpr const double MZ    = 91.1876;
	constexpr const double PI = std::numbers::pi;
	constexpr const double PI_2 = PI*PI;
	constexpr const double PI_3 = PI*PI*PI;
	constexpr const double Zeta2 = PI*PI/6.0;
	constexpr const double Zeta3 = 1.2020569031595942854;
	/** @} */

	// TODO: better handle defaults
	/**
	 *  @defgroup defaults Program Defaults
	 *  @{
	 */
	constexpr const uint DISTS = 37;
	constexpr const uint INTERP_POINTS = 4;
	constexpr const uint DEFAULT_ITERATIONS = 10;
	constexpr const uint DEFAULT_TRUNC_IDX = 5;
	/** @{ */

	/**
	 *  @brief Template class for more easily typing an @a std::vector with multiple layers of nesting.
	 */
	template <typename T, uint N>
	struct MultiDimVector
	{
		typedef typename MultiDimVector<T,N-1>::type Nested;
		typedef std::vector<Nested> type;
	};
	template <typename T>
	struct MultiDimVector<T,1>
	{
		typedef std::vector<T> type;
	};

	/** Enum for defining a set of standard logging types. */
	enum LogType : int
	{
		LOG_DEBUG = 0,
		LOG_INFO,
		LOG_WARNING,
		LOG_ERROR,
		LOG_ERROR_NOQUIT,
		LOG_NUM_LOG_TYPES
	};
	inline std::array<std::string_view, LOG_NUM_LOG_TYPES> log_string_reps{"DEBUG", "INFO", "WARNING", "ERROR", "ERROR"};

	// TODO: better handle setting global flags
	inline bool debug_flag{false}; //!< whether to print logged messaged tagged with LOG_DEBUG
	/** Getter/Setter for the debug flag */
	inline bool& getDebugFlag() { return debug_flag; }
	
	inline bool use_log_output_stream{false}; //!< flag for whether to use an additional logging output stream
	inline std::reference_wrapper<std::ostream> log_output_stream = std::ref(std::cout); //!< additional output stream
	/** Setter for the additional logging output stream */
	inline void set_log_output_stream(std::ostream& os) { log_output_stream = os; use_log_output_stream = true; }

	/**
	 *  @brief Prints a message to standard out and possibly an additional stream with a nice prefix.
	 *  @param log_type a type as in the enum @a LogType
	 *  @param prefix an identifier of some kind to include within the message, often a function name
	 *  @param fmt_string the format string of the log message, as used in std::format
	 *  @param args args used to format the string, as used in std::format
	 */
	template <typename... TArgs>
	void log(uint log_type, std::string_view prefix, std::format_string<TArgs...> fmt_string, TArgs&& ...args)
	{
		if (log_type == LOG_DEBUG && !debug_flag)
			return;
			
		std::string log_text = std::vformat(fmt_string.get(), std::make_format_args(args...));
		std::string all_text = std::format("[{}] {}: {}\n", log_string_reps[log_type], prefix, log_text);
		if (use_log_output_stream)
			log_output_stream.get() << all_text;
		std::cout << all_text;

		if (log_type == LOG_ERROR)
			exit(EXIT_FAILURE);
	}

	/**
	 *  @brief a concept to require that a type has a @a value_type as well as begin and end iterators
	 */
	template <typename TContainer>
	concept CContainer = requires(TContainer&& t)
	{
		typename TContainer::value_type;
	    { std::begin(t) } -> std::input_or_output_iterator;
		{ std::end(t) } -> std::sentinel_for<decltype(std::begin(t))>;
	};

	/**
	 *  @brief Returns a string with the values of the container separated by a comma and a space
	 *  @param vec The vector to turn into a string
	 */
	template <CContainer TContainer>
	std::string vec_to_str(TContainer const& vec)
	{
		using value_type = TContainer::value_type;
		std::ostringstream ss{};
		std::copy(vec.begin(), vec.end(), std::ostream_iterator<value_type>(ss, ", "));
		return std::move(ss.str());
	}
}; // namespace Candia2


#endif // __COMMON_HPP
