/*---------------------------------------------------------------------------*/
/*Version 26-March-1997                                      File: sunalpha.c*/
/*---------------------------------------------------------------------------*/
/* Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica             */
/* Universidad Complutense de Madrid, 28040-Madrid, Spain                    */
/* E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es                    */
/*---------------------------------------------------------------------------*/
/* This program is free software; you can redistribute it and/or modify it   */
/* under the terms of the GNU General Public License as published by the Free*/
/* Software Foundation; either version 2 of the License, or (at your option) */
/* any later version. See the file gnu-public-license.txt for details.       */
/*---------------------------------------------------------------------------*/
/*Comment                                                                    */
/*                                                                           */
/*Program: sunalpha                                                          */
/*Classification: input/output                                               */
/*Description: Transforms REDUCEME images from sun architecture to alpha     */
/*architecture (and viceversa).                                              */
/*                                                                           */
/*Comment                                                                    */
/*---------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>

void invert4(FILE *fp1, FILE *fp2, unsigned int n);
void direct1(FILE *fp1, FILE *fp2, unsigned int n);
unsigned int invert4read(FILE *fp1, FILE *fp2, unsigned int opc);

int main(int argc, char *argv[])
{
  FILE *fp1,*fp2;
  unsigned int opc;
  unsigned int nchar;
  unsigned char byte1,byte2,byte3,byte4;

  printf("*************************************************");
  printf("******************************\n");
  printf("REDUCEMEv4.0                  Welcome to sunalpha");
  printf("        Version: 27-March-1997\n");
  printf("-------------------------------------------------");
  printf("------------------------------\n\n");

  if(argc != 3) 
    {
    printf("ERROR: remeber usage --> sunalpha file_input file_output\n");
    exit(1);
    }
  printf("(1) file_input  corresponds to alpha architecture\n");
  printf("    file_output corresponds to sun architecture\n\n");
  printf("(2) file_input  corresponds to sun architecture\n");
  printf("    file_output corresponds to alpha architecture\n\n");
  do 
    {
    printf("Option (1/2)? ");
    scanf("%d",&opc);
    if ((opc != 1) && (opc != 2))
      printf("ERROR: invalid entry. Try again.\n");
    } 
  while ((opc != 1) && (opc != 2));

  if( (fp2 = fopen(argv[2],"rb")) != NULL)
    {
    printf("ERROR: file_output %s already exist\n",argv[2]);
    exit(1);
    }

  if( (fp1 = fopen(argv[1],"rb")) == NULL)
    {
    printf("ERROR: while opening %s\n",argv[1]);
    exit(1);
    }

  if( (fp2 = fopen(argv[2],"wb")) == NULL)
    {
    printf("ERROR: while opening %s\n",argv[2]);
    exit(1);
    }

  /*CLAVE*/
  invert4(fp1,fp2,1);  /* numero de bytes de CLAVE: 12 = 0000000C */
  direct1(fp1,fp2,12); /* cadena */
  invert4(fp1,fp2,1);  /* numero de bytes de CLAVE: 12 = 0000000C */

  /*NSCAN,NCHAN*/
  invert4(fp1,fp2,4); /* numero de bytes de NSCAN,NCHAN: 8 = 00000008 */
                      /* NSCAN,NCHAN: 2 x 4 bytes */
                      /* numero de bytes de NSCAN,NCHAN: 8 = 00000008 */

  /*STWV,DISP*/
  invert4(fp1,fp2,4); /* numero de bytes de STWV,DISP: 8 = 00000008 */ 
                      /* STWV,DISP: 2 x 4 bytes */
                      /* numero de bytes de STWV,DISP: 8 = 00000008 */

  /*AIRMASS*/
  invert4(fp1,fp2,3); /* numero de bytes de AIRMASS: 4 = 00000004 */ 
                      /* AIRMASS: 1 x 4 bytes */
                      /* numero de bytes de AIRMASS: 4 = 00000004 */

  /*TIMEXPOS*/
  invert4(fp1,fp2,3); /* numero de bytes de TIMEXPOS: 4 = 00000004 */ 
                      /* TIMEXPOS: 1 x 4 bytes */
                      /* numero de bytes de TIMEXPOS: 4 = 00000004 */

  /*TrueLen of OBJECT*/
  invert4(fp1,fp2,1); /* numero de bytes de TrueLen of OBJECT: 4 = 00000004 */
  nchar=invert4read(fp1,fp2,opc); /* TrueLen of OBJECT */
  invert4(fp1,fp2,1); /* numero de bytes de TrueLen of OBJECT: 4 = 00000004 */
  printf("TrueLen of OBJECT..: %d\n",nchar);
  /*OBJECT*/
  invert4(fp1,fp2,1); /* numero de bytes de OBJECT = TrueLen of OBJECT */
  direct1(fp1,fp2,nchar); /* OBJECT: nchar bytes */
  invert4(fp1,fp2,1); /* numero de bytes de OBJECT = TrueLen of OBJECT */

  /*TrueLen of FITSFILE*/
  invert4(fp1,fp2,1); 
  nchar=invert4read(fp1,fp2,opc); 
  invert4(fp1,fp2,1);
  printf("TrueLen of FITSFILE: %d\n",nchar);
  /*FITSFILE*/
  invert4(fp1,fp2,1);
  direct1(fp1,fp2,nchar);
  invert4(fp1,fp2,1);

  /*TrueLen of COMMENT*/
  invert4(fp1,fp2,1);
  nchar=invert4read(fp1,fp2,opc); 
  invert4(fp1,fp2,1);
  printf("TrueLen of COMMENT.: %d\n",nchar);
  /*OBJECT*/
  invert4(fp1,fp2,1);
  direct1(fp1,fp2,nchar);
  invert4(fp1,fp2,1);

  do 
    {
    byte1=getc(fp1);
    byte2=getc(fp1);
    byte3=getc(fp1);
    byte4=getc(fp1);
    if (!feof(fp1)) putc(byte4,fp2);
    if (!feof(fp1)) putc(byte3,fp2);
    if (!feof(fp1)) putc(byte2,fp2);
    if (!feof(fp1)) putc(byte1,fp2);
    } 
  while (!feof(fp1));

  
  if( fclose(fp1) != 0) 
    {
    printf("ERROR: while closing %s\n",argv[1]);
    exit(1);
    }  

  if( fclose(fp2) != 0)
    {
    printf("ERROR: while closing %s\n",argv[2]);
    exit(1);
    }  
  return(0);
}

/*****************************************************************************/
/* invierte el orden de 4 bytes n veces y escribe en el fichero de salida */
void invert4(FILE *fp1, FILE *fp2, unsigned int n)
{
  unsigned int i;
  unsigned char byte1,byte2,byte3,byte4;

  for (i=1; i<=n; i++) {
    byte1=getc(fp1);
    byte2=getc(fp1);
    byte3=getc(fp1);
    byte4=getc(fp1);
    putc(byte4,fp2);
    putc(byte3,fp2);
    putc(byte2,fp2);
    putc(byte1,fp2);
    }
}

/*****************************************************************************/
/* escribe en el fichero de salida n bytes */
void direct1(FILE *fp1, FILE *fp2, unsigned int n)
{
  unsigned int i;
  unsigned char byte;

  for (i=1; i<=n; i++) {
    byte=getc(fp1);
    putc(byte,fp2);
    }
}

/*****************************************************************************/
/* invierte el orden de 4 bytes, escribe en el fichero de salida y devuelve */
/* el entero correspondiente */
unsigned int invert4read(FILE *fp1, FILE *fp2, unsigned int opc)
{
  unsigned int num;
  unsigned char byte1,byte2,byte3,byte4;

  byte1=getc(fp1);
  byte2=getc(fp1);
  byte3=getc(fp1);
  byte4=getc(fp1);
  putc(byte4,fp2);
  putc(byte3,fp2);
  putc(byte2,fp2);
  putc(byte1,fp2);
  if(opc == 1) /* alpha */
    {
    num= byte1+byte2*256+byte3*65536+byte4*16777216;
    }
  else         /* sun */
    {
    num= byte4+byte3*256+byte2*65536+byte1*16777216;
    }
  return(num);
}
