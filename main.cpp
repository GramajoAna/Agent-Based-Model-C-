/* *************************************************************************************************************
			Mosquito program for one block, Fortran version (2021). Author Fabiana Laguna.
			Extended program for a grid L x L,  C++ version (2021). Author Ana A. Gramajo.
							/\
							||
We developed an agent-based model (ABM) composed of a collection of autonomous decision-making individuals following a set of rules. Agents are embedded in a two-dimensional grid where we set a spatial distribution of water. 
//*************************************************************************************************************
*/

#include <cstdlib>
#include <iostream>
#include <fstream>		
#include <vector>		
#include <functional>   
#include <algorithm>    
#include <cmath>

#include "cpu_timer.h"
#include "mosquitas.h"


int main(){

//************************  open files ********************************************************************
	std::ofstream outfile,outfile1,outfile2;
    	outfile.open("Total-population.dat");
    	outfile1.open("Adult-population.dat");
    	outfile2.open("Aquatic-population.dat");

//***************** CPU timer opened ************************************************************************   	
	cpu_timer Reloj_CPU;
	Reloj_CPU.tic();

//***************** Initializing struct bichos **************************************************************
	bichos mosquitas(Ninicial);

//***************** Iterate through dia (days = 1; days < 400) **********************************************
	for(int dia = 1; dia <= Ndias; dia ++)
	{
	
		int tovip=tiempo_entre_oviposiciones(dia);
		mosquitas.mortalidades(dia);
			if(dia > 120 && dia < 320)mosquitas.descacharrado(dia);
		mosquitas.conteo_huevos(dia);
		mosquitas.reproducir(dia,tovip);
//************* Calculate the total population, adult population and aquatic population *********************
		int vivas=mosquitas.vivos(dia);
		int adultos=mosquitas.adultos(dia);
		int acuaticos=mosquitas.acuaticos(dia);
//**************************** outputs files *****************************************************************
		outfile << dia << "\t" << vivas << std::endl;
		outfile1 << dia << "\t" << adultos << std::endl;
		outfile2 << dia << "\t" << acuaticos << std::endl;
//***********************************************************************************************************
		mosquitas.envejecer(dia); //depende del orden
		mosquitas.delay(dia); //depende del orden
	}

//**************** CPU timer closed ([t]= min)*************************************************************************   	
    double t=Reloj_CPU.tac()/60000; 
    printf("Tiempo en CPU: %lf minutos\n",t);

//close files
outfile.close();
outfile1.close();
outfile2.close();

	return 0;							

}
