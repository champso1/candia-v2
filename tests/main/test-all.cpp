#include "Candia-v2/Common.hpp"
using namespace Candia2;

#include <vector>
#include <stdexcept>
#include <filesystem>
#include <format>
#include <string>
#include <string_view>
#include <ranges>
#include <algorithm>
#include <iterator>
namespace fs = std::filesystem;

#include "util.hpp"

static std::vector<fs::path> create_datafiles(
	std::vector<int> orders,
	std::vector<std::string_view> names,
	int iterations, int trunc_idx,
	std::vector<double> mur2_muf2_vals
);
static fs::path find_compare_file(fs::path const& candia_file);
static bool compare_files(fs::path const& candia_file, fs::path const& other_file);
static void validate_file(fs::path const& candia_file);

int main()
{
	auto orders = std::vector{0, 1, 2, 3};
	auto names = std::vector{"lo", "nlo", "nnlo", "n3lo"};
	auto iterations = 15;
	auto trunc_idx = 10;
	auto mur2_muf2_vals = std::vector{0.5, 1.0, 2.0};

	auto datafile_list = create_datafiles(
		{0, 1, 2, 3},
	    {"lo", "nlo", "nnlo", "n3lo"},
		13, 10,
		{0.5, 1.0, 2.0});

	std::ranges::for_each(datafile_list, [](auto&& datafile){ validate_file(datafile); });
}



static std::vector<fs::path> create_datafiles(
	std::vector<int> orders,
	std::vector<std::string_view> names,
	int iterations, int trunc_idx,
	std::vector<double> mur2_muf2_vals
) {
	std::vector<fs::path> datafile_list{};
	
	for (auto order : orders) {
		auto name = names[order];
		for (auto mur2_muf2 : mur2_muf2_vals) {
			auto mur2_muf2_str =
				mur2_muf2 == 0.5 ? "05" :
				mur2_muf2 == 1.0 ? "10" :
				mur2_muf2 == 2.0 ? "20" : throw std::runtime_error("invalid mur2_muf2 value");
			for (auto use_trunc : {true, false}) {
				auto use_trunc_str = use_trunc ? "1" : "0";

				auto datafile_name = std::format(
					"{}-{}-{}",
					name,
					use_trunc ? "trunc" : "exact",
					std::format("mu{}", mur2_muf2_str)
				);

				auto datafile_path = fs::path("./data")/fs::path(datafile_name).replace_extension(".dat");
				datafile_list.emplace_back(datafile_path);
				if (fs::exists(datafile_path)) {
					log(LOG_INFO, "test-all.cpp", "datafile '{}' already exists. skipping...", datafile_path.string());
					continue;
				} else {
					log(LOG_INFO, "test-all.cpp", "datafile '{}' does not exist. creating...", datafile_path.string());
				}

				auto command = std::format(
					"./evolve {} {} {} {} {} 0 {}",
					order,
					iterations,
					trunc_idx,
					mur2_muf2,
					use_trunc_str,
					datafile_name
				);

				log(LOG_INFO, "test-all.cpp", "{}", command);
				std::system(command.c_str());
			}
		}
	}

	assert(!datafile_list.empty(), "test-all.cpp: returned 0 datafiles");
	return datafile_list;
}


static fs::path find_compare_file(fs::path const& candia_path)
{
	using std::operator""sv;
	static constexpr auto delim{"-"sv};

	auto get_token_as_sv = [](auto view, auto pos) {
		auto it = view.begin();
		while (pos-->0)
			it++;
		return std::string_view((*it).begin(), (*it).end());
	};
	
	auto candia_datafile_filename = candia_path.filename().replace_extension().string();
	auto view = std::views::split(candia_datafile_filename, delim);
	
	auto order_name = get_token_as_sv(view, 0);
	auto mur2_muf2_str = get_token_as_sv(view, 2);

	if (
		(
			(order_name.compare("nnlo"sv) == 0) || (order_name.compare("n3lo"sv) == 0)
		) && (
			mur2_muf2_str.compare("mu10"sv) != 0
		)
	) {
		return fs::path{};
	} else if (mur2_muf2_str.compare("mu10"sv) == 0){
		auto p = fs::path("./hoppet-data")/fs::path(std::format("out-{}.dat", order_name));
		if (!fs::exists(p)) {
			log(LOG_ERROR, "test-all.cpp",
				"selected '{}' (to compare with '{}'), but it doesn't exist",
				p.string(), candia_path.string());
		}
		return p;
	} else {
	    auto p = fs::path("./pegasus-data")/fs::path(std::format("out-{}-{}.dat", order_name, mur2_muf2_str));
		if (!fs::exists(p)) {
			log(LOG_ERROR, "test-all.cpp",
				"selected '{}' (to compare with '{}'), but it doesn't exist",
				p.string(), candia_path.string());
		}
		return p;
	}

	
}

static bool compare_files(fs::path const& candia_file, fs::path const& other_file)
{
	auto candia_file_res = read_candia_file(candia_file, 37);
	auto [other_xtab, other_dists_raw]  = read_other_file(other_file, 37);

	if (!std::ranges::equal(candia_file_res.xtab, other_xtab)) {
		log(LOG_ERROR_NOQUIT, "test-all.cpp",
			"  - '{}' vs '{}': mismatch in xtab ({} vs {})",
			candia_file.string(), other_file.string(), candia_file_res.xtab.size(), other_xtab.size()
		);
		return false;
	}

	auto candia_dists_fixed = fix_dists(candia_file_res.dists_ntabbed, 0);
	auto other_dists_fixed  = fix_dists(other_dists_raw, 0);

	auto xtab = std::span(candia_file_res.xtab.begin()+1, candia_file_res.xtab.end()-1);
	auto candia_dists = [&](auto j) {
		auto const& arr_full = candia_dists_fixed[j];
		return std::span(arr_full.begin()+1, arr_full.end()-1);
	};
	auto other_dists = [&](auto j) {
		auto const& arr_full = other_dists_fixed[j];
		return std::span(arr_full.begin()+1, arr_full.end()-1);
	};

	auto diff_func = [](auto candia, auto other) {
		auto avg = (candia+other)/2.0;
		return std::abs((candia-other)/avg);
	};

	auto dists_view =
		std::views::iota(0, 13) |
		std::views::filter(
			[](auto i) {
				if (i==6 || i==12)
					return false;
				return true;
			});
	auto dists = std::vector(dists_view.begin(), dists_view.end());
	auto dists_enumerate =
		std::views::iota(uint{0}, dists.size()) |
		std::views::transform([&](auto i){ return std::make_pair(i, dists[i]); });

	assert(
		(
			xtab.size() == candia_dists(0).size()
		) && (
			candia_dists(0).size() == other_dists(0).size()
		) && (
			candia_dists_fixed.size() == other_dists_fixed.size()
		) && (
			candia_dists_fixed.size() == dists.size()
		),
		"test-all.cpp", "failed to setup restricted range for dists");
	
	auto xtab_enumerate =
		std::views::iota(uint{0}, xtab.size()) |
		std::views::transform([&](auto i){ return std::make_pair(i, xtab[i]); });

	auto failures = std::vector<std::pair<uint,uint>>{};

	static constexpr auto cutoff = 1.0e-4; // equivalent to saying percent < 0.01%
	for (auto [ji, j] : dists_enumerate) {
		for (auto [xi,x] : xtab_enumerate) {
			auto candia_val = candia_dists(ji)[xi];
			auto other_val = other_dists(ji)[xi];
			auto diff = diff_func(candia_val, other_val);
			if (!(diff < cutoff))
				failures.emplace_back(std::make_pair(ji,xi));
		}
	}

	if (failures.empty())
		return true;

	log(LOG_ERROR_NOQUIT, "test-all.cpp",
		"Between '{}' and '{}:",
		candia_file.string(), other_file.string());
	for (auto [ji, xi] : failures) {
		auto j = dists[ji];
		auto x = xtab[xi];
		auto candia_val = candia_dists(ji)[xi];
		auto other_val = other_dists(ji)[xi];
		auto diff = diff_func(candia_val, other_val);

		log(LOG_ERROR_NOQUIT, "test-all.cpp",
			"  - PDF {} does not agree at x={} ({} vs {}), diff={:.4e}",
			j, x, candia_val, other_val, diff);
	}
	
	return false;
}

static void validate_file(fs::path const& candia_file)
{
	auto other_file = find_compare_file(candia_file);
	if (other_file.empty()) {
		log(LOG_INFO, "test-all.cpp", "file '{}' didn't have a valid comparison file. skipping...", candia_file.string());
		return;
	} else {
		log(LOG_INFO, "test-all.cpp", "choosing '{}' to compare with '{}'",
			other_file.string(), candia_file.string());
	}

	auto res = compare_files(candia_file, other_file);
	if (res)
		log(LOG_INFO, "test-all.cpp", "'{}' is good!", candia_file.string());
	else
		log(LOG_INFO, "test-all.cpp", "'{}' failed. See above", candia_file.string());
}
