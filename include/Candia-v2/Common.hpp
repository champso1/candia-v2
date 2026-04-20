/**
 *  @file Common.hpp
 *  @brief Contains type definitions, constants, logging, and other misc definitions.
 */

#pragma once

#include <functional>
#include <string_view>
#include <unordered_map>
#include <vector>
#include <numbers>
#include <format>
#include <iostream>
#include <sstream>
#include <iterator>
#include <ranges>
#include <mutex>
#include <string_view>
#include <array>
#include <filesystem>
#include <fstream>

using uint = unsigned;

#define NUM_SINGLET_THREADS 8

namespace Candia2
{
	/**
	*  @defgroup constants Constants
	*  @{
	*/
	constexpr double CF    = 4.0/3.0;
	constexpr double NC    = 3.0;
	constexpr double TR    = 0.5;
	constexpr double MZ    = 91.1876;
	constexpr double PI = std::numbers::pi;
	constexpr double PI_2 = PI*PI;
	constexpr double PI_3 = PI*PI*PI;
	constexpr double Zeta2 = PI*PI/6.0;
	constexpr double Zeta3 = 1.2020569031595942854;
	/** @} */

	/**
	*  @defgroup defaults Program Defaults
	*  @{
	*/
	constexpr const uint DISTS = 37;
	constexpr const uint INTERP_POINTS = 4;
	constexpr const uint DEFAULT_ITERATIONS = 10;
	constexpr const uint DEFAULT_TRUNC_IDX = 5;
	constexpr const uint NUM_SUBTRACT_PDFS = 2;
	/** @{ */


	// colors for printing to the terminal
	inline constexpr char const* ANSI_COLOR_RED =         "\x1b[31m";
	inline constexpr char const* ANSI_COLOR_GREEN =       "\x1b[32m";
	inline constexpr char const* ANSI_COLOR_YELLOW =      "\x1b[33m";
	inline constexpr char const* ANSI_COLOR_BLUE =        "\x1b[34m";
	inline constexpr char const* ANSI_COLOR_MAGENTA =     "\x1b[35m";
	inline constexpr char const* ANSI_COLOR_CYAN =        "\x1b[36m";
	inline constexpr char const* ANSI_COLOR_RESET =       "\x1b[0m";
	inline constexpr char const* ANSI_LINEFEED_CLEAR =    "\033[2K";
	inline constexpr auto ANSI_LINEFEED_UP   = [](uint count){ return std::format("\033[{}F", count); };
	inline constexpr auto ANSI_LINEFEED_DOWN = [](uint count){ return std::format("\033[{}E", count); };
	inline constexpr std::string_view loading_block("█");

	
	/** Enum for defining a set of standard logging types. */
	enum LogType : int
	{
		LOG_DEBUG = 0,
		LOG_INFO,
		LOG_WARNING,
		LOG_ERROR,
		LOG_ERROR_NOQUIT,
		LOG_THREAD,
		LOG_NUM_LOG_TYPES
	};
	inline std::array<std::string_view, LOG_NUM_LOG_TYPES> log_string_reps{
		"DEBUG",
		"INFO",
		"WARNING",
		"ERROR",
		"ERROR",
		"THREAD"};
	inline std::array<std::string_view, LOG_NUM_LOG_TYPES> log_string_colors{
		ANSI_COLOR_GREEN,
		ANSI_COLOR_RESET,
		ANSI_COLOR_YELLOW,
		ANSI_COLOR_RED,
		ANSI_COLOR_RED,
		ANSI_COLOR_CYAN};

	/** @brief struct to store flags/options for logging */
	struct LogOptions final
	{
		bool silent{false}; //!< suppresses all messages
		bool show_debug_messages{false}; //!< switch for showing debug messages. only useful for debugging, default off
		bool show_thread_output{false};  //!< switch for showing the output from threads. often floods the console, default off

		bool use_log_output_stream{false}; //!< switch for whether we are logging to another output stream, like a file
		std::reference_wrapper<std::ostream> log_output_stream{std::ref(std::cerr)}; //!< actual output stream. only used if @a use_log_output_stream is true

		inline static LogOptions makeDefault()
		{
			return LogOptions{};
		}
	};
	/** @brief global options struct */
	inline LogOptions _log_options = LogOptions::makeDefault();
	/** @brief returns the global options for setting values */
	inline LogOptions& getLogOptions() { return _log_options; }
	/** @brief setter for apply some desired options */
	inline void setLogOptions(LogOptions const& o) { _log_options = LogOptions{o}; };

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
		if (getLogOptions().silent)
			return;
		if (log_type == LOG_DEBUG && !getLogOptions().show_debug_messages)
			return;
		if (log_type == LOG_THREAD && !getLogOptions().show_thread_output)
			return;
			
		std::string log_text = std::vformat(fmt_string.get(), std::make_format_args(args...));
		std::string all_text = std::format("{}[{}] {}: {}{}\n",
			log_string_colors[log_type], log_string_reps[log_type], prefix, log_text, ANSI_COLOR_RESET);
		if (getLogOptions().use_log_output_stream)
			getLogOptions().log_output_stream.get() << all_text;
		std::cout << all_text;

		if (log_type == LOG_ERROR)
			exit(EXIT_FAILURE);
	}

	/**
	 *  @brief A simple wrapper for printing to standard out and possibly an additional stream.
	 *  @param fmt_string the format string of the log message, as used in std::format
	 *  @param args args used to format the string, as used in std::format
	 */
	template <typename... TArgs>
	void log(std::format_string<TArgs...> fmt_string, TArgs&&... args)
	{
		std::string log_text = std::vformat(fmt_string.get(), std::make_format_args(args...));
		if (getLogOptions().use_log_output_stream)
			getLogOptions().log_output_stream.get() << log_text;
		std::cout << log_text;
	}

	inline std::unordered_map<uint, uint> log_threads_line_offset{};
	inline void registerThreadLogs(std::vector<uint> const& ids)
	{
		for (uint i=0; i<ids.size(); ++i) {
			log_threads_line_offset[ids[i]] = ids.size()-i;
			std::cout << '\n';
		}
	}
	inline void unregisterThreadLogs([[maybe_unused]] std::vector<uint> const& ids)
	{
		log_threads_line_offset.clear();
	}
	inline std::mutex log_threads_mutex;
	/**
	 *  @brief Prints messages to standard out for threaded iterations in a nice way that doesn't flood stdout
	 *  @param thread_idx used in threaded printing, 
	 *  @param val the i in i/j, for the iteration count
	 *  @param end the j in i/j, for the iteration count
	 *  @param prefix an identifier of some kind to include within the message, often a function name
	 */
	template <typename... TArgs>
	void logThreadIterations(uint thread_idx, uint val, uint end, std::string_view prefix)
	{
		if (!getLogOptions().show_thread_output)
			return;
		auto log_type = LOG_THREAD;

		auto count = log_threads_line_offset[thread_idx];
		auto ansi_jump_up   = ANSI_LINEFEED_UP(count);
		auto ansi_jump_down = ANSI_LINEFEED_DOWN(count);
		double ratio = static_cast<double>(val)/static_cast<double>(end);
		int num_blocks = static_cast<int>(ratio*50.0);
		std::string blocks{};
		for (int i=0; i<num_blocks; ++i)
			blocks += std::string(loading_block);

		std::string all_text = std::format("{}[{}] {}: [{}] Iteration {:0>2}/{} ({: >3}%) [{: <50}]",
			log_string_colors[log_type],
			log_string_reps[log_type], prefix,
			thread_idx, val, end, static_cast<uint>(ratio*100.0),
			blocks,
			ANSI_COLOR_RESET);
		if (getLogOptions().use_log_output_stream)
			getLogOptions().log_output_stream.get() << all_text;

		{
			std::lock_guard<std::mutex> guard{log_threads_mutex};

			std::cout << ansi_jump_up << ANSI_LINEFEED_CLEAR
					  << all_text << std::flush
					  << ansi_jump_down << std::flush;
		}
	}

	inline void startLogIterations(){ return; }
	inline void endLogIterations(){ std::cout << '\n'; }
	/**
	 *  @brief Prints messages to standard out for iterations in a nice way that doesn't flood stdout
	 *  @param val the i in i/j, for the iteration count
	 *  @param end the j in i/j, for the iteration count
	 *  @param prefix an identifier of some kind to include within the message, often a function name
	 */
	template <typename... TArgs>
	void logIterations(uint val, uint end, std::string_view prefix)
	{
		if (!getLogOptions().show_thread_output)
			return;
		auto log_type = LOG_INFO;

		double ratio = static_cast<double>(val)/static_cast<double>(end);
		int num_blocks = static_cast<int>(ratio*50.0);
		std::string blocks{};
		for (int i=0; i<num_blocks; ++i)
			blocks += loading_block;

		std::string all_text = std::format("\r{}[{}] {}: Iteration {:0>2}/{} ({: >3}%) [{: <50}]{}",
			log_string_colors[log_type],
			log_string_reps[log_type], prefix,
			val, end, static_cast<uint>(ratio*100.0),
			blocks,
			ANSI_COLOR_RESET);
		if (getLogOptions().use_log_output_stream)
			getLogOptions().log_output_stream.get() << all_text;
		std::cout << all_text << std::flush;
	}

	template <typename... TArgs>
	void assert(bool stmnt, std::format_string<TArgs...> fmt_string, TArgs&&... args)
	{
		if (!stmnt)
			log(LOG_ERROR, "ASSERT", fmt_string, std::forward<TArgs>(args)...);
	}

	/**
	 *  @brief a concept to require that a type is a range or a view, i.e. a generic container of values
	 */
	template <typename TContainer>
	concept RangeContainer = std::ranges::range<TContainer> || std::ranges::view<TContainer>;

	/**
	 *  @brief Returns a string with the values of the container separated by a comma and a space
	 *  @param vec The vector to turn into a string
	 */
	template <RangeContainer TContainer>
	std::string vec_to_str(TContainer const& vec, std::string const& delim=", ")
	{
		using value_type = decltype(*std::ranges::begin(vec));
		std::ostringstream ss{};
		std::ranges::copy(vec, std::ostream_iterator<value_type>(ss, delim.c_str()));
		return ss.str();
	}

	/**
	 *  @brief Returns a string with the values of the container separated by a comma and a space
	 *  @param vec The vector to turn into a string
	 *  this version, opposed to the @a vec_to_str, will not leave a trailing delimiter
	 */
	template <RangeContainer TContainer>
	std::string vec_to_str2(TContainer const& vec, std::string const& delim=", ")
	{
	    std::ostringstream ss{};
		auto it = std::ranges::begin(vec);
		auto end = std::ranges::end(vec);

		ss << std::setprecision(10) << std::scientific;
		if (it != end) {
			ss << *it;
			++it;
		}

		for (; it != end; ++it)
			ss << delim << *it;

		return ss.str();
	}

	extern thread_local int thread_index;
	inline void initializeThreadIndex(int index)
	{
		thread_index = index;
	}

	/**
	 *  @brief returns a string of the file contents
	 */
	inline std::string read_file(std::filesystem::path const& filepath)
	{
		if (!std::filesystem::exists(filepath))
			log(LOG_ERROR, "read_file()", "Failed to open file '{}'", filepath.string());

		std::ifstream infile(filepath);
		return std::string(
			std::istreambuf_iterator<char>(infile),
			std::istreambuf_iterator<char>{});
	}
}; // namespace Candia2
