// Old conditions for the singlet stuff in the Array/Function Grid
// evolution/resummation in Evolve2
// ====================================================================================================
/*
if (resum_threshold)
{
    for (uint j=0; j<=1; ++j)
        temp_arr_singlet[j*31] = _S2[0][j][0];
}
else if (resum_tab)
{
    for (uint j=0; j<=1; ++j)
        _F2[j*31] = _S2[0][j][0];
}
*/
/*
if (resum_tab)
{
    _S2[0][1][1][k] = RecRelS_1(_S2[0][1][0], k, getSplitFunc("P0qq")) + RecRelS_1(_S2[0][0][0], k, getSplitFunc("P0qg"));
    _S2[0][0][1][k] = RecRelS_1(_S2[0][1][0], k, getSplitFunc("P0gq")) + RecRelS_1(_S2[0][0][0], k, getSplitFunc("P0gg"));

    for (uint t=0; t<=_trunc_idx; ++t)
    {
        for (uint j=0; j<=1; j++)
            _F2[j*31][k] += _S2[t][j][1][k] * std::pow(_alpha1, t) * std::pow(L1, n)/Factorial(n);
    }
}
else if (resum_threshold)
{
    _S2[0][1][1][k] = RecRelS_1(_S2[0][1][0], k, getSplitFunc("P0qq")) + RecRelS_1(_S2[0][0][0], k, getSplitFunc("P0qg"));
    _S2[0][0][1][k] = RecRelS_1(_S2[0][1][0], k, getSplitFunc("P0gq")) + RecRelS_1(_S2[0][0][0], k, getSplitFunc("P0gg"));

    for (uint t=0; t<=_trunc_idx; ++t)
    {
        for (uint j=0; j<=1; j++)
            temp_arr_singlet[j][k] += _S2[t][j][1][k] * std::pow(_alpha1, t) * std::pow(L1, n)/Factorial(n);
    }
}
*/
/* 
for (uint n=1; n<_iterations; ++n)
{
    std::cerr << "[DGLAP2] Singlet iteration " << n << '\n';
    for (uint k=0; k<_grid.Size()-1; ++k)
    {
        _S2[0][1][n][k] = RecRelS_1(_S2[0][1][n-1], k, getSplitFunc("P0qq")) + RecRelS_1(_S2[0][0][n-1], k, getSplitFunc("P0qg"));
        _S2[0][0][n][k] = RecRelS_1(_S2[0][1][n-1], k, getSplitFunc("P0gq")) + RecRelS_1(_S2[0][0][n-1], k, getSplitFunc("P0gg"));
    }	
}
std::cerr << "[DGLAP2] Finished singlet evolution\n";

std::cerr << "[DGLAP2] Performing singlet resummation.\n";
if (resum_tab)
{
    for (uint j=0; j<=1; j++)
    {
        _F2[j*31] = _S2[0][j][0];
        std::cerr << "[DGLAP2] Singlet distribution " << j << '\n';
        for (uint k=0; k<_grid.Size()-1; k++)
        {
            for (uint n=1; n<_iterations; n++)
            {
                for (uint t=0; t<=_trunc_idx; t++)
                    _F2[j*31][k] += _S2[t][j][n][k] * std::pow(_alpha1, t)*std::pow(L1, n)/Factorial(n);
            }
        }
    }
}
else if (resum_threshold)
{
    for (uint j=0; j<=1; j++)
    {
        for (uint k=0; k<_grid.Size()-1; k++)
        {
            for (uint n=1; n<_iterations; n++)
            {
                for (uint t=0; t<=_trunc_idx; t++)
                    _S2[0][j][0][k] += _S2[t][j][n][k]*std::pow(_alpha1, t)*std::pow(L1, n)/Factorial(n);
            }
        }
    }
}
std::cerr << "[DGLAP2] Finished singlet resummation.\n";
*/
// ====================================================================================================
// Old conditions for the non-singlet stuff in the Array/Function Grid
// evolution/resummation in Evolve2
/*
if (resum_threshold)
{
    for (uint j=13; j<=30+_nf; ++j)
    {
        temp_arr[j] = _A2[j][0];
        if (j == (12+_nf))
            j = 31;
        
    }
}
else if (resum_tab)
{
    for (uint j=13; j<=30+_nf; ++j)
    {
        _F2[j] = _A2[j][0];
        if (j == (12+_nf))
            j = 31;
    }
}
*/
/*
if (resum_tab)
{
    for (uint j=13; j<=12+_nf; j++)
    {
        _A2[j][1][k] = RecRelLO(_A2[j][0], k, getSplitFunc("P0ns"));
        _F2[j][k] += _A2[j][1][k]*std::pow(L1, n)/Factorial(n);
    }

    for (uint j=32; j<=30+_nf; j++)
    {
        _A2[j][1][k] = RecRelLO(_A2[j][0], k, getSplitFunc("P0ns"));
        _F2[j][k] += _A2[j][1][k]*std::pow(L1, n)/Factorial(n);
    }
    
}
else if (resum_threshold)
{
    for (uint j=13; j<=12+_nf; j++)
    {
        _A2[j][1][k] = RecRelLO(_A2[j][0], k, getSplitFunc("P0ns"));
        temp_arr[j][k] += _A2[j][1][k]*std::pow(L1, n)/Factorial(n);
    }

    for (uint j=32; j<=30+_nf; j++)
    {
        _A2[j][1][k] = RecRelLO(_A2[j][0], k, getSplitFunc("P0ns"));
        temp_arr[j][k] += _A2[j][1][k]*std::pow(L1, n)/Factorial(n);
    }
}
*/