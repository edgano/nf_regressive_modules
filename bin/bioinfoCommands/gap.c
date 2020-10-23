#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>

/*
    Compile: gcc -Wall gap.c -o gap
    gap <input_file> <boolean>
                        true: fasta -> gap
                        false: gap -> fasta
*/
char* inputFasta; //= "test.fa";
char* outputFileFA2Gap = "resultF2G.gap";
char* outputFileGap2FA = "resultG2F.fa";
FILE *outFileFa2Gap;
FILE *outFileGap2Fa;

bool gaps;
char* fasta2Gap;
char gap = '-';

int writeChar2File(FILE *outFile, char char2write){
    int i;
    if (outFile != NULL) {
        //printf ("printChar: %c \n", char2write);
        fputc (char2write, outFile);		
        //printf("done printChar\n");
   }
   else {
  	  printf("\n Unable to Create or Open the %s File",outFile);
   }
   return 0;
}
void printString(FILE *filePoint, char *string2print){
    int len = strlen(string2print);
    int i;
    //printf ("printString: %s \n", string2print);
    for (i = 0; i < len; i++) {
        writeChar2File(filePoint,string2print[i]);
	}
}
void convert_final_digit(int digit){
    int i;
    for( i = 0; i < digit; i = i + 1 ){
        writeChar2File(outFileGap2Fa,gap);
    }
}
void convertGaps (char *line){
    char *c;
    char aux[100];  //number of digits for gaps... posibility of ERROR
    char nn; 
    int numberGap=0;
    int lineLenght=strlen(line);
    int lineIndex=1;

    for (c=line; *c; c++) {
        nn = *c;
        if (nn == gap)  {            // gap found
            numberGap++;
            if(lineLenght==lineIndex){   //if it's last postion of line -> printString
                sprintf(aux, "%d", numberGap);  //convert int to char*
                //printf("... ....%s\n", aux);
                printString(outFileFa2Gap,aux);
                //reset counter
                numberGap=0;
            }
        }else{            
            if (numberGap!=0){      // 1st non-gap after gap
                //print number of gaps
                sprintf(aux, "%d", numberGap);  //convert int to char*
                //printf("... ....%s\n", aux);
                printString(outFileFa2Gap,aux);
                //reset counter
                numberGap=0;
            }
            //print nucleotide
            //printf("... NN ... %c\n", nn);
            writeChar2File(outFileFa2Gap,nn);
        }
        lineIndex++;
    }
}
void convertFasta (char *line){
    char *c;
    char nn; 
    int check=0;
    int auxDigit;
    int lineLenght=strlen(line);
    int lineIndex=1;

    for (c=line; *c; c++) {
        nn = *c;
        if (!isdigit(nn) || check==1){  //it's a nucleotide or non the first digit
            if(isdigit(nn)){
                //printf("Digists in a row >> %c\n",nn);  
                //second or more digits in a row
                //nothing to do... it's already captured in the ATOI()
                //printf("DIGIT> lenght %d - index %d\n",lineLenght,lineIndex);

                if(lineLenght==lineIndex){   //if it's last postion of line -> printString
                    convert_final_digit(auxDigit);
                }
            }else{
                if (check==0){  //digit already covnerted
                //printf("nucleotide: %c\n",nn);
                //****printf("%c",nn);
                writeChar2File(outFileGap2Fa,nn);
                }else{          //first time to nucelotide after digit
                    //printf("auxDigit: : : %i",auxDigit);
                    convert_final_digit(auxDigit);
                    check=0;
                    //printf("total gaps: %d \nfirst neuclotide %c\n", auxDigit,nn);
                    //****printf("%c",nn);
                    writeChar2File(outFileGap2Fa,nn);
                }
            }
        }else{              //its the FIRST digit
            //it graps all the digirs in teh auxDigit 
            //but the pointer in the char* not... it will
            //go in the if bellow with CHECK==1...
            auxDigit= atoi(c);
            //printf("<< digit %d\n",auxDigit);
            check=1;
            //TODO -> check if its just 0 --> ERROR
            //printf("ELSE>lenght %d - index %d\n",lineLenght,lineIndex);
            if(lineLenght==lineIndex){   //if it's last postion of line -> printString
                convert_final_digit(auxDigit);
            }
        }   
        lineIndex++;
    }
}
void readFasta(char *fileName){
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
 
    fp = fopen(fileName, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);
 
    int state = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        /* Delete trailing newline */
        if (line[read - 1] == '\n')
            line[read - 1] = 0;
        /* Handle comment lines*/
        if (line[0] == '>') {
            if (state == 1){
                //printf("\n");
            }
            //****printf(">%s: \n", line+1);
            if (gaps){  //choose wich is the output file
                printString(outFileFa2Gap,">");
                printString(outFileFa2Gap,line+1);
                printString(outFileFa2Gap,"\n");
            }else{
                printString(outFileGap2Fa,">");
                printString(outFileGap2Fa,line+1);
                printString(outFileGap2Fa,"\n");
            }
            state = 1;
        } else {
            /* Print everything else */
            if (gaps){  //fasta >> gaps
                convertGaps(line);
                printString(outFileFa2Gap,"\n");
            }else {     // gaps >> fasta
                //printf("convert2Fasta: %s \n", line);
                convertFasta(line);
                printString(outFileGap2Fa,"\n");
            }
        }
    }
    fclose(fp);
    if (line)
        free(line);
}
void convert2Fasta(char *filename){
    gaps = false;
    readFasta(filename);
}
void convert2Gap(char *filename){
    gaps = true;
    readFasta(filename);
}
FILE *openResultFile(char *fileName){
    FILE *fileAddress = fopen(fileName, "w");
    return fileAddress;
}
int closeResultFile(FILE *outFile){
    // printf("\nWe have written the file successfully\n");
	fclose(outFile);	
}
void strip_ext(char *fname){
    char *end = fname + strlen(fname);

    while (end > fname && *end != '.' && *end != '\\' && *end != '/') {
        --end;
    }
    if ((end > fname && *end == '.') && (*(end - 1) != '\\' && *(end - 1) != '/')) {
        *end = '\0';
    }  
}
int main(int argc, char* argv []){

    inputFasta=argv[1];             // file to convert
    fasta2Gap=argv[2];            //second arg ==true if we want fasta -> gap format

    if (strcmp(fasta2Gap, "true") == 0){     // FASTA to GAP
        //open result files
        outFileFa2Gap = openResultFile(outputFileFA2Gap);
        printf("## FASTA TO GAP ##\n    Input file: %s\n",inputFasta);

        convert2Gap(inputFasta);
        closeResultFile(outFileFa2Gap);

        printf("    Output file: %s\n",outputFileFA2Gap);

    }else if(strcmp(fasta2Gap, "false") == 0){              //GAP to FASTA
        //open result files
        outFileGap2Fa = openResultFile(outputFileGap2FA);
        printf("## GAP TO FASTA ##\n    Input file: %s\n",inputFasta);

        convert2Fasta(inputFasta);
        closeResultFile(outFileGap2Fa);

        printf("    Output file: %s\n",outputFileGap2FA);
    }else {
        printf("Not a valid mode to parse %s\nUse \"true\" to convert Fa->Gap and \"false\" to go Gap->Fa\n",fasta2Gap);
    }
    return 0;
}