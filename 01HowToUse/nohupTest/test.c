#include<stdio.h>
#include <unistd.h>



int main() {
   /* my first program in C */
   printf("Hello, World! \n");
	for(int i = 1;i<1000;i++)
	{
	//	printf("i = %d\n",i);
		// sleep(1);
	//	wait(1);
		printf("i = %d HELLO\n",i);
		fflush(stdout);
		sleep(3);
	}   
   return 0;
}
