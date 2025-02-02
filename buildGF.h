void buildGF(int *gf, short int n, int polinom)//
{
    for ( int i=1;i<n;i++)
    {
        gf[i]=gf[i-1]>>1;
        if((gf[i-1]%2)==1)
            {
                gf[i]^=polinom;
            }
    }
}