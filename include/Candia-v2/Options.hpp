#ifndef __OPTIONS_HPP
#define __OPTIONS_HPP

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

#endif // __OPTIONS_HPP
