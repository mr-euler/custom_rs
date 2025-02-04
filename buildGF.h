void buildGF(int *gf,int *reverse_gf, int gf_size, int polinom)//
{
    for ( int i=1;i<gf_size;i++)
    {
        gf[i]=gf[i-1]>>1;
        if((gf[i-1]%2)==1)
            {
                gf[i]^=polinom;
            }
		reverse_gf[gf[i]]=i;
    }
    reverse_gf[gf[0]]=0;
	reverse_gf[0]=-1;

}
