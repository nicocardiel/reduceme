/*---------------------------------------------------------------------------*/
/*Version 07-September-2007                         File: interchange_32_64.c*/
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
/*Program: interchange_32_64                                                 */
/*Classification: input/output                                               */
/*Description: Transforms REDUCEME images from 32 bits to 64 bits            */
/*architecture (and viceversa).                                              */
/*                                                                           */
/*Comment                                                                    */
/*---------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>

void separator32(FILE *fp_input, FILE *fp_output);
void separator64(FILE *fp_input, FILE *fp_output);

int main(int argc, char *argv[])
{
  //declaracion de variables
  unsigned int i,icheck;
  unsigned int nbits;
  FILE *fp_input,*fp_output;
  unsigned char byte;

  if(argc != 3)
  {
    printf("*************************************************");
    printf("******************************\n");
    printf("REDUCEMEv4.0          Welcome to interchange_32_64");
    printf("   Version: 07-September-2007\n");
    printf("-------------------------------------------------");
    printf("------------------------------\n");
    printf("Usage: interchange_32_64 input_filename output_filename\n");
    exit(1);
  }

  //---------------------------------------------------------------------------
  //comprobamos que es un fichero en formato REDUCEME
  //primero probamos con una arquitectura de 32 bits
  if((fp_input=fopen(argv[1],"rb")) == NULL)
  {
    printf("ERROR while opening the input file \"%s\".\n",argv[1]);
    exit(1);
  }
  for(i=1;i<=4;i++)
    byte=getc(fp_input); //saltamos 4 primeros bytes
  icheck=0;
  for(i=1;i<=12;i++)
  {
    if(icheck==0)
    {
      byte=getc(fp_input);
      if((unsigned int)byte != 96+i) icheck=1;
    }
  }
  if(icheck != 0) //no es REDUCEME de 32 bits; probamos con 64 bits
  {
    fclose(fp_input);
    fp_input=fopen(argv[1],"rb");
    for(i=1;i<=8;i++)
      byte=getc(fp_input); //saltamos 8 primeros bytes
    icheck=0;
    for(i=1;i<=12;i++)
    {
      if(icheck==0)
      {
        byte=getc(fp_input);
        if((unsigned int)byte != 96+i) icheck=1;
      }
    }
    if(icheck !=0) //no es REDUCEME de 64 bits; fin del programa
    {
      printf("ERROR: the file \"%s\" does not seem to have REDUCEME format.\n",argv[1]);
      exit(1);
    }
    else //es REDUCEME de 64 bits
    {
      nbits=64;
    }
  }
  else //es REDUCEME de 32 bits
  {
    nbits=32;
  }
  fclose(fp_input);
  fp_input=fopen(argv[1],"rb");

  //---------------------------------------------------------------------------
  //abrimos el fichero de salida
  if((fp_output=fopen(argv[2],"rb")) != NULL)
  {
    printf("ERROR: the output file \"%s\" already exists.\n",argv[2]);
    exit(1);
  }
  if((fp_output=fopen(argv[2],"wb")) == NULL)
  {
    printf("ERROR while opening the output file \"%s\".\n",argv[2]);
    exit(1);
  }

  //---------------------------------------------------------------------------
  //transformamos fichero segun el caso
  if(nbits==32) //pasamos de 32 a 64 bits
  {
    do
    {
      separator32(fp_input,fp_output);
    }
    while(!feof(fp_input));
  }
  else if(nbits==64) //pasamos de 64 a 32 bits
  {
    do
    {
      separator64(fp_input,fp_output);
    }
    while(!feof(fp_input));
  }
  else
  {
    printf("ERROR: nbits=%d.\n",nbits);
    exit(1);
  }

  //---------------------------------------------------------------------------
  //cerramos los ficheros
  if(fclose(fp_input) != 0)
  {
    printf("ERROR while closing input file.\n");
    exit(1);
  }
  if(fclose(fp_output) != 0)
  {
    printf("ERROR while closing output file.\n");
    exit(1);
  }

  //fin del programa
  return(0);
}

/*****************************************************************************/
/* Lee un bloque formado por SEPARADOR32+DATOS+SEPARADOR32, donde SEPARADOR32
 * esta formado por 4 bytes y nos indica el numero de bytes que constituyen
 * DATOS. Escribe SEPARADOR64+DATOS+SEPARADOR64, donde SEPARADOR64 esta
 * formado por 8 bytes.*/
void separator32(FILE *fp_input, FILE *fp_output)
{
  unsigned char byte0,byte1,byte2,byte3,byte4;
  unsigned long i,ndata;

  byte0=0;

  //leemos separador de 4 bytes
  byte1=getc(fp_input);
  byte2=getc(fp_input);
  byte3=getc(fp_input);
  byte4=getc(fp_input);
  //escribimos un separador de 8 bytes
  if(!feof(fp_input))
  {
    putc(byte1,fp_output);
    putc(byte2,fp_output);
    putc(byte3,fp_output);
    putc(byte4,fp_output);
    putc(byte0,fp_output);
    putc(byte0,fp_output);
    putc(byte0,fp_output);
    putc(byte0,fp_output);

    //calculamos numero de bytes en los datos
    ndata=byte1+byte2*256+byte3*65536+byte4*16777216;

    //leemos los datos y los escribimos
    if(ndata > 0)
    {
      for (i=1;i<=ndata;i++)
      {
        byte1=getc(fp_input);
        if(!feof(fp_input)) putc(byte1,fp_output);
      }
    }

    //leemos separador de 4 bytes
    byte1=getc(fp_input);
    byte2=getc(fp_input);
    byte3=getc(fp_input);
    byte4=getc(fp_input);
    //escribimos un separador de 8 bytes
    if(!feof(fp_input))
    {
      putc(byte1,fp_output);
      putc(byte2,fp_output);
      putc(byte3,fp_output);
      putc(byte4,fp_output);
      putc(byte0,fp_output);
      putc(byte0,fp_output);
      putc(byte0,fp_output);
      putc(byte0,fp_output);
    }
  }
}


/*****************************************************************************/
/* Lee un bloque formado por SEPARADOR64+DATOS+SEPARADOR64, donde SEPARADOR64
 * esta formado por 8 bytes y nos indica el numero de bytes que constituyen
 * DATOS. Escribe SEPARADOR32+DATOS+SEPARADOR32, donde SEPARADOR32 esta
 * formado por 4 bytes.*/
void separator64(FILE *fp_input, FILE *fp_output)
{
  unsigned char byte1,byte2,byte3,byte4,byte5,byte6,byte7,byte8;
  unsigned long i,ndata;

  //leemos separador de 8 bytes
  byte1=getc(fp_input);
  byte2=getc(fp_input);
  byte3=getc(fp_input);
  byte4=getc(fp_input);
  byte5=getc(fp_input);
  byte6=getc(fp_input);
  byte7=getc(fp_input);
  byte8=getc(fp_input);
  if(!feof(fp_input))
  {
    if((byte5!=0)||(byte6!=0)||(byte7!=0)||(byte8!=0))
    {
      printf("ERROR: byte5, byte6, byte7 or byte8 != 0.\n");
      exit(1);
    }
  }
  //escribimos un separador de 4 bytes
  if(!feof(fp_input))
  {
    putc(byte1,fp_output);
    putc(byte2,fp_output);
    putc(byte3,fp_output);
    putc(byte4,fp_output);

    //calculamos numero de bytes en los datos
    ndata=byte1+byte2*256+byte3*65536+byte4*16777216;

    //leemos los datos y los escribimos
    if(ndata > 0)
    {
      for (i=1;i<=ndata;i++)
      {
        byte1=getc(fp_input);
        if(!feof(fp_input)) putc(byte1,fp_output);
      }
    }

    //leemos separador de 8 bytes
    byte1=getc(fp_input);
    byte2=getc(fp_input);
    byte3=getc(fp_input);
    byte4=getc(fp_input);
    byte5=getc(fp_input);
    byte6=getc(fp_input);
    byte7=getc(fp_input);
    byte8=getc(fp_input);
    if(!feof(fp_input))
    {
      if((byte5!=0)||(byte6!=0)||(byte7!=0)||(byte8!=0))
      {
        printf("ERROR: byte5, byte6, byte7 or byte8 != 0.\n");
        exit(1);
      }
    }
    //escribimos un separador de 4 bytes
    if(!feof(fp_input))
    {
      putc(byte1,fp_output);
      putc(byte2,fp_output);
      putc(byte3,fp_output);
      putc(byte4,fp_output);
    }
  }
}
