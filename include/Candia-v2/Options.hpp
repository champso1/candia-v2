/**
 *  @file Options.hpp
 *  @brief Contains the @a Options class which handles the description of options for a class in a consistent manner
 */

#pragma once

#include <concepts>

namespace Candia2
{
    template <typename TOptions>
    requires std::default_initializable<TOptions>
    struct OptionsBase 
    {
        using options_type = TOptions;
        options_type options{};
        virtual options_type& getOptions() { return options; };
    };
}
