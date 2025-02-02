int binaryToDecimal(int num,  int *p) 
{ 
    int count=0;
    int dec_value = 0; 
    int base = 1;
    int last_digit = num % 10; 
    while (num) { 
      //3   int last_digit = num % 10; 
        num /= 10; 
        dec_value += last_digit * base; 
        base *= 2;
        last_digit = num % 10;
        count++;

    } 
    
    *p=count;
    return dec_value; 
}