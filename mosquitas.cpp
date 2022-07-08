/* *************************************************************************************************************
			Mosquito program for one block, Fortran version (2021). Author Fabiana Laguna.
			Extended program for a grid L x L,  C++ version (2021). Author Ana A. Gramajo.
							/\
							||
We developed an agent-based model (ABM) composed of a collection of autonomous decision-making individuals following a set of rules. Agents are embedded in a two-dimensional grid where we set a spatial distribution of water. 
//*************************************************************************************************************


The simulation setup:
	- Spatiality. We consider a square grid of size L x L = M, the number of blocks of the simulated city. 
	- Initial conditions.  We define the total number of mosquitoes in the system N(t=0), i.e., initially we set:
		State[i]: 0 (alive).
		Age[i]: a random number between 19 and 25 days (extracted from a uniform distribution).
    		Bucket[i]: i. 
    		Block[i]: a random number between 0 and M (extracted from a uniform distribution).
    		Adulthood[i]: the pupa matures between 15 and 19 days (extracted from a uniform distribution).
    		LifeSpan[i]: a random number between 27 and 32 days (extracted from a uniform distribution).

We initially consider one mosquito per bucket. 
We calculate the daily mosquito population along a year. For simplicity, we consider 400 days for the simulation. 

Daily we computed the adult and aquatic population considering:
	- Oviposition. Following the Bergero2012 results, we assume that the female lays between 10 and 35 eggs per oviposition (extracted from a uniform distribution). 
	- For the number of ovipositions that each female can have, we consider the model given by Otero2006: 
		- one at 18°C,
		- four or five at 23°C and 27°C, 
		- six at 30°C.
	- Mortality rates. We consider the daily mortality of eggs, pupae, larvae, and adults independent of temperature and density. The data was extracted from Otero2006 where:
        	- The mortality of the eggs is chosen to be 0.01 1/day.
        	- The death of the larvae is approximated by 0.01 1/day.
        	- The intrinsic mortality of a pupa has been considered as 0.01 1/day.
        	- The daily mortality in the pupal stage associated with the unsuccessful emergence of the adult individual is assume a mortality of 0.171/day.
We consider that the mortality of adults is 0.01 1/day according to the number of days that each individual lives.
	- Death by old age. The mosquito dies when it reaches its chosen life span.
	- Seasons. We discretize the data corresponding to the mean daily temperatures for the period July 2001–July 2002 used by Otero2008 from Buenos Aires, Argentina. The simulation starts on July 1st and, the temperature at different seasons is defined as follows:
    		-  T=18°C in the range of days [1,80) and (320,400]. 
    		-  T=23°C between [80, 140] days. 
    		-  T=30°C within (140,260) days.
    		-  T=27°C between  [260,320] days. 
	- Saturation of the buckets. We suppose a maximum number of eggs allowed per bucket (800). 
	- Bucket transfer. In the present model, we assume that adult females can move with a higher probability to another bucket into the same block (80 %) and only can move with a lower probability to the nearest neighbor (20 %).
	- Advertising campaigns. We introduce the discarding effectiveness through the emptying of a percentage (PROP) of the buckets per block during the hottest months, i.e., December, January, and February. In our model, this season corresponds to the range of days (120,320). 
	- Frequency of discarding. We chose from a uniform distribution in the interval [1,14] for each container (unless otherwise indicated). In this way, all the buckets ARE SYNCHRONIZED To emptying at once.
	- Time-delay in the bucket disponibility (nTau). Once the buckets are emptied during the advertising campaigns, we add a time-delay measured in days until being available.
	- Aquatics hibernation. We include a delay on the development of eggs, larval and pupal stages in the coolest days, i.e., between [1,80) and (320,400) days. 
 	- Daily count. We count daily the number of adult and aquatic individuals.
	- Daily mosquito population aging. We increase the age of the mosquitoes by one day.   
*/

#include <cstdlib>
#include <iostream>
#include <fstream>		
#include <vector>		
#include <functional>   
#include <algorithm>    
#include <cmath>

#include "cpu_timer.h"
#include "ran2.h"
#include "parameters.h" //for input parameters

int tiempo_entre_oviposiciones(int dia)
	{
	int t;
		if(dia >= 1 && dia < 80) t = tovip3;      
		if(dia >= 80 && dia <= 140)t = tovip2b;   
  		if(dia > 140 || dia < 260) t = tovip1;    
 		if(dia >= 260 && dia <= 320)t = tovip2a;  
 		if(dia > 320) t = tovip3;                 

	return t;
	}

struct bichos{

	std::vector<int> estado;  
	std::vector<int> edad;    
	std::vector<int> tacho;   
	std::vector<int> TdV; 
	std::vector<int> pupacion; 
	std::vector<int> manzana; 
	std::vector<int> tach; 
	std::vector<int> Tau; 

	std::vector<int> N_mobil; 

	std::vector<float> E; 
	std::vector<int> manzana_del_tacho; 
	std::vector<std::vector<int> > tachos_por_manzana; //NUEVO tachos_por_manzana[i]=vector de tachos de manzana i 
	std::vector<int> NroTachos; //NUEVO
	std::vector<float> N_descach; // NUEVO array donde almaceno los tachos que se van a eliminar en cada manzana


	//construct	
	bichos(int N_)
	{

	// alocamos el maximo posible
	estado.resize(MAXIMONUMEROBICHOS);	
	tacho.resize(MAXIMONUMEROBICHOS);	
	edad.resize(MAXIMONUMEROBICHOS);
	pupacion.resize(MAXIMONUMEROBICHOS);
	TdV.resize(MAXIMONUMEROBICHOS);
	manzana.resize(MAXIMONUMEROBICHOS);
	tach.resize(MAXIMONUMEROBICHOS);
	Tau.resize(MAXIMONUMEROBICHOS);

	N_mobil.resize(1);

	E.resize(NUMEROMANZANAS); 
	tachos_por_manzana.resize(NUMEROMANZANAS); 
	manzana_del_tacho.resize(MAXIMONUMEROBICHOS);
	NroTachos.resize(NUMEROMANZANAS);
	N_descach.resize(NUMEROMANZANAS);

	std::fill(estado.begin(),estado.end(),0);
	std::fill(edad.begin(),edad.end(),0);
	std::fill(tacho.begin(),tacho.end(),0);
	std::fill(pupacion.begin(),pupacion.end(),0);
	std::fill(TdV.begin(),TdV.end(),0);
	std::fill(manzana.begin(),manzana.end(),0);
	std::fill(Tau.begin(),Tau.end(),0);

	std::fill(E.begin(),E.end(),prop);	


	
	//****************************  Initial conditions ***************************************************************************
 	for(int i=0;i < N_;i++)
 	{
 		tacho[i] = i;			                               //tacho en el que se encuentra la mosquita
			if(DISTRIBUCIONTACHOS == 1)
			{
				manzana_del_tacho[tacho[i]] = int (i/5);
			}		
			else{
			manzana_del_tacho[tacho[i]] = int(ran2(&semilla)*NUMEROMANZANAS);
			}  
			
		manzana[i] = manzana_del_tacho[tacho[i]];                //manzana en la que está el tacho i
		tachos_por_manzana[manzana[i]].push_back(tacho[i]);      //para identificar los tachos tengo en la manzana    	
		int j = tacho[i];
		Tau[j] = 0;   						   //NUEVOS, todos los tachos disponibles
		int tachosxmanzana = tachos_por_manzana[manzana[i]].size();//nro de tachos x manzana
		
			if(tachosxmanzana <= 9)
			{ 
				estado[i] = ESTADOVIVO; 	      		             		
				edad[i] = ran2(&semilla)*7+19; 	              
				pupacion[i] = tpupad-2+(ran2(&semilla)*5);            
				TdV[i] = ran2(&semilla)*6+27;	                     
			}
	}

//***************** verificacion del llenado de tachos x manz. y de la disponibilidad de tacho x manz.************************

	int nmanzanas=tachos_por_manzana.size();

		for(int i=0;i < nmanzanas;i ++){
			int ntachos = tachos_por_manzana[i].size();   //nro de tachos por manzana
			int contTach = 0;                             //contador para el nro de tachos por manzana

			for(int j = 0;j < ntachos;j ++)
			{
				//std::cout << (tachos_por_manzana[i])[j] << "\n";
				contTach ++;
				NroTachos[i] = contTach;  
			}//cierro for para ntachos
			//std::cout << "\nNro de tachos en la manzana\t" << contTach << "\n";
		}//cierro for para nmanzanas
	
       //almaceno en el array N_descach el número de tachos que se van a descacharrar en cada manzana
        std::transform(NroTachos.begin(), NroTachos.end(), E.begin(), N_descach.begin(),std::multiplies<float>());
      
	std::cout << "Inicializacion" << std::endl;

	N_mobil[0]=N_;
	};//close constructor	

//*************************************************** Mortalities calculation ***************************************************	
	void mortalidades(int dia)
	{
	int N=N_mobil[0];
		for(int i=0;i < N;i++)
		{
			// aquatics mortality = eggs + larvae + pupae mortalitites
			if (estado[i] == ESTADOVIVO && edad[i] < pupacion[i]){if(ran2(&semilla) < moracu) estado[i] = ESTADOMUERTO;}
			// emergence
			if (estado[i] == ESTADOVIVO && edad[i] == pupacion[i]){if(ran2(&semilla) < morpupad) estado[i] = ESTADOMUERTO;}
			// adult mortality
			if (estado[i] == ESTADOVIVO && edad[i] > pupacion[i]){if(ran2(&semilla) < morad) estado[i] = ESTADOMUERTO;} 
			// the oldest mosquitoes died 
			if (estado[i] == ESTADOVIVO && edad[i] >= TdV[i]) estado[i] = ESTADOMUERTO;
		} 
	};
	
//*************************************************** Buckets discarding ***************************************************	
	void descacharrado(int dia){
	int N=N_mobil[0];

	int nmanzanas=tachos_por_manzana.size();	
	int azar,nTau;

	//Frecuency discarding input in paramaters.h
	if(FRECUENCIA==1)
	{
		azar = 7; // 7 days between two successive discards
	}
	else
	{
		azar = 1 + ran2(&semilla)*14;// a uniform distribution of discarding frequencies between 1 and 14 days
	}
	
	//The effectiveness of advertising campaigns by emptying a given number of buckets per block during the hot and temperate seasons
	//in main() corresponds to 120 < dia < 320
		if(dia%azar == 0)
		{
   			for(int i = 0;i < nmanzanas; i ++)	//loop over the number of blocks
   			{
   		        int NroDescach = round(N_descach[i]);             //number of buckets to dicard in the block i
   		        int ntachos = tachos_por_manzana[i].size();       //number of buckets in the block i

   				for(int itach = 0; itach < NroDescach; itach ++) //loop over the number of buckets to discard in the block i
   				{
   			        	int n = ran2(&semilla)*ntachos;    	//sort the index for the bucket to discard
   		 	        	int ntach = (tachos_por_manzana[i])[n];//bucket selected to discard 
   		 	        	
   		 	        		//conditional for the time-delay in the bucket disponibility input from parameters.h
						if(DELAYTACHO == 1)
						{
							nTau = 10; //fixed
						}	
						else
						{
							nTau=1+ran2(&semilla)*20;//a uniform distribution of time-delay between 1 and 20 days
						}	
					Tau[ntach]=nTau;// change the disponibility for the bucket discarding

					//change the state to death for all the mosquitoes that belong to the bucket discarding selected
						for(int i=0;i < N;i++)
						{
		  					if (estado[i] == ESTADOVIVO && edad[i] < pupacion[i] && tacho[i] == ntach)estado[i]=ESTADOMUERTO;
						}
				}    	
   			}
		}	   
	};

	//count the number of acuatic population in each buckets
	void conteo_huevos(int dia)
	{
	int N=N_mobil[0];
		for(int i = 0;i < N; i++)
		{ 
			tach[i] = 0;
   			if(edad[i] < pupacion[i] && estado[i] == ESTADOVIVO){ 
    			int j = tacho[i]; 
	    		tach[j] ++;
			}
		} 
	};


	//We calculate the (new) block considering their nearest neighbor blocks where the female oviposit.
	
	int sorteo_manzana_vecina(int manz)
	{

		std::vector<int> manzanas_vecinas(5); //we considerer an arrays dimension where the index corresponds to the nearest neighbor blocks of the block input (manz) and the input (manz)
				
		int x = manz%LADO;
		int y = int(manz/LADO);
		manzanas_vecinas[0] = (x-1+LADO)%LADO+LADO*y; 	//left 
		manzanas_vecinas[1] = (x+1+LADO)%LADO+LADO*y; 	//right
		manzanas_vecinas[2] = LADO*((y-1+LADO)%LADO)+x; 	//down 
		manzanas_vecinas[3] = LADO*((y+1+LADO)%LADO)+x; 	//up
		manzanas_vecinas[4] = manz; 				//center
		
		int nvecinos = 6; 					//counting the center an their neighbors
		int r = int((rand()*1.0/RAND_MAX)*nvecinos); 		//randon number between 0 and 5
		int manzana_sorteada;					//for the new block

		// conditional: if r == 1 the mosquito stays in their same block (center), if r > 1 then moving to another block with index equals to r		
		if(r > 1)
		{
			manzana_sorteada = manz;
		}
		else
		{		
			manzana_sorteada = manzanas_vecinas[r];
		}

		return manzana_sorteada;
	};

//*************************************************** Bucket transfer ***************************************************	
	void reproducir(int dia,int tovip)
	{	

	int indice=N_mobil[0]; //number of mosquitoes
	
	int estamanzana=0;

		for(int i=0;i < indice;i++) //loop over the number of mosquitoes
		{
			//conditional if the mosquito is alive and adult to oviposit 
			if(estado[i] == ESTADOVIVO && edad[i] > pupacion[i] && edad[i]%tovip == 0)
			{
	
/*Buckets transfer local (mosquitoes can move to another bucket within the same block) and global (mosquitoes can also move to another bucket of the first neighbouring block)*/
				//check if the bucket where live the mosquito is saturated
				if (tach[tacho[i]] < sat) estamanzana = manzana_del_tacho[tacho[i]];//block of the bucket[i]
				
				int manzanadeltacho = sorteo_manzana_vecina(estamanzana);//call a function that returns the block where the mosquito oviposit (could be in the same block or nearby)
				
				int ntachosNew = tachos_por_manzana[manzanadeltacho].size();//number of buckets in the (new) block 
				int nNew = ran2(&semilla)*ntachosNew;//sort a new number between 0 and number of buckets in the (new) block 
				int k=(tachos_por_manzana[manzanadeltacho])[nNew];	      //the new bucket

				//check if the new bucket is saturated and available 
					if (tach[k] < sat && Tau[k] == 0)
					{
				    		int iovip = 10 + (ran2(&semilla)*25); //sort a number of eggs that female oviposit in the new bucket
   					    		for(int ik = 0; ik < iovip; ik ++)// for each new mosquito set the initial conditions
   					    		{ 
 								estado[indice] = ESTADOVIVO;// state alive
 								edad[indice] = 1;   //age = 1 day
 								tacho[indice] = k; //bucket
		         					pupacion[indice] = tpupad-2+(ran2(&semilla)*5);//adulthood	
	         						TdV[indice] = ran2(&semilla)*6+27;  //lifespan
								manzana[indice] = manzana_del_tacho[tacho[indice]];//block
								int j = tacho[indice];
								Tau[j] = 0; //bucket disponibility
		 						tach[j] ++; //counting the new ones in the bucket
								indice ++; //add the new ones 	
   							} 
					}
					//if the bucket is saturated then looking for other located at the (new) block
					else 
					{
						int l=0;
					//while over in all the buckets	
							while(l < ntachosNew) 
							{ 
								int t=(tachos_por_manzana[manzanadeltacho])[l]; //l is the index for the new bucket
							//we check if the new bucket is saturated and available and repeat the conditional line 334
			         	 				if(t!=k && tach[t] < sat && Tau[t] == 0){
						    				int iovip = 10 + (ran2(&semilla)*25); 
   							    			for(int ik = 0; ik < iovip; ik ++)
   							    			{ 
 							        			estado[indice] = ESTADOVIVO;
 							        			edad[indice] = 1;   
 							        			tacho[indice] = t; 
		         				    				pupacion[indice] = tpupad-2+(ran2(&semilla)*5);	
				         						TdV[indice] = ran2(&semilla)*6+27;  
											manzana[indice] = manzana_del_tacho[tacho[indice]];
											int j = tacho[indice];
											Tau[j] = 0; 
 											tach[j] ++;
											indice ++;
   										}
							    			l = ntachosNew; //flag
					    				}
								l ++;
							}	
					} 				
   			} 
		}

		//Refresing the number of mosquitoes 
		N_mobil[0]=indice;

	};
		
	// We calculate the total population
	int vivos(int dia){

	int N = N_mobil[0];
	int poblacion = 0;

  		for(int i = 0; i < N; i++)
  		{
			if (estado[i] == ESTADOVIVO)poblacion++;
		}	
	
	return poblacion;
	};

	//We calculate the Aquatic population
	int acuaticos(int dia){

	int N = N_mobil[0];
  	int ac = 0;  
  		for(int i = 0; i < N; i++)
  		{
			if (estado[i] == ESTADOVIVO && edad[i] < pupacion[i])ac++; 
		}

	return ac;
	};

	//We calculate the adult population
	int adultos(int dia){

	int N = N_mobil[0];
  	int ad = 0;  
  	
	  	for(int i = 0; i < N; i++){
			if (estado[i] == ESTADOVIVO && edad[i] >= pupacion[i])ad++; //adultos tot	
		}

	return ad;

	};

	//We increased the mosquitos age daily
	void envejecer(int dia){
	int N = N_mobil[0];
	
  		for(int i = 0; i < N; i++){
  		
  			//In winter, the aquatic population doesn't old
			if(dia < 80 || dia > 320)
			{
    				if(estado[i] == ESTADOVIVO && edad[i] > pupacion[i]) edad[i]++;
    			} 
			else
			{
				if(estado[i] == ESTADOVIVO)edad[i]++;
			} 

		}
	};

	//we reduce in 1 day the bucket disponibility
	void delay(int dia){
	
	int N = N_mobil[0];

  		for(int i = 0; i < N; i++){ 
			if(Tau[i] > 0) Tau[i]--;
		}
	};

};// close struct bichos



int main(){

//************************ Input files ********************************************************************
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
//**************************** outputs file *****************************************************************
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
