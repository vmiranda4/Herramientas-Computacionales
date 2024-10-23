#include <stdio.h>
#include <stdlib.h>

int main(){

// reading to a file
FILE * vFile;
vFile = fopen("var.dat","r");
int j;
double x=0;
double y[3];
for(j=1; j<4; j++){
fscanf(vFile, "%lf", &x);
y[j-1]=x;
}
fclose(vFile);


}