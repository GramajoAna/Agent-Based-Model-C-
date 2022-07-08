

//-----------------------------------------------------------------------------------------------------------------------
//					Initial conditions
// Ninicial : total number of mosquitoes in the system
// NUMEROTACHOS : we initially consider one mosquito per bucket so Ninicial == NUMEROTACHOS
// MAXIMONUMEROBICHOS : max arrays dimension
// ESTADOMUERTO : dead mosquito
// ESTADOVIVO: live mosquito
// DISTRIBUCIONTACHOS: number of buckets per block. 
//-----------------------------------------------------------------------------------------------------------------------

#define Ninicial	    	125	
#define NUMEROTACHOS		125	
#define MAXIMONUMEROBICHOS	9000000 
#define DISTRIBUCIONTACHOS	1  // Uniform distribution: 1 (five buckets per block), Random distribution: 0 (until 9 buckets per block)

//State of the mosquito
#define ESTADOMUERTO		1	
#define ESTADOVIVO		0	


//Number of days days for the simulation
#define Ndias			400


//-----------------------------------------------------------------------------------------------------------------------
//					Spaciallity
//-----------------------------------------------------------------------------------------------------------------------
// GRID L x L 
#define LADO           	5       	
#define NUMEROMANZANAS		LADO*LADO   	//number of blocks 

//-----------------------------------------------------------------------------------------------------------------------
//					Social Parameters
// prop: discarding effectiveness (between 0 and 1)
// FRECUENCIA: frecuency of dicarding
// DELAYTACHO: Delay in the recovery of the buckets
//-----------------------------------------------------------------------------------------------------------------------

#define prop 			0.6	
#define FRECUENCIA		1	//1 para frecuencia semanal de descacharrado, 0 para frecuencia de descacharrado cada 7 días en promedio
#define DELAYTACHO		1	//1 para delay de los tachos de 10 días, 0 para delay de 10 días en promedio


//-----------------------------------------------------------------------------------------------------------------------
//					Biologial Parameters
// morhue: eggs mortality (Otero 2006)
// morlar: larva mortality (Otero 2006)
// morpup: death of pupae (Otero 2006)
// morad: adults mortality 
// moracu: aquatics mortality = eggs + larva + pupae mortalities
// morpupad: emergence 
// Oviposition. According to the model presented by Otero 2006, the average number of eggs laid by one adult female in one oviposition is 63, while the field data obtained by Bergero 2012 gives a number of eggs much lower. Following the Bergero et. al results, we assume that the female lays between 10 and 35 eggs per oviposition (extracted from a uniform distribution). Moreover, for the number of ovipositions that each female can have, we consider the model given by Otero 2006: one at 18◦C, four or five at 23◦C and 27◦C, and six at 30◦C.
//-----------------------------------------------------------------------------------------------------------------------

// Mortalities [1/days] 
#define morhue 		0.01    
#define morlar 		0.01     
#define morpup 		0.01    
#define morad 			0.01 	
#define moracu 		0.03	
#define morpupad 		0.17	

// Times [days] 

#define tpupad	 		17  	// when the pupae became adult  
#define tovip1 		2	// time between two ovipositions at T = 27 °C
#define tovip2a 		3  	// time between two oviposiciones at T = 23 °C
#define tovip2b		4   	// time between two oviposiciones at T = 23 °C
#define tovip3 		30	// time between two oviposiciones at T = 18 °C
#define sat 			800     // maximum number of eggs allowed per bucket

// Seed for random generator ran2() 
long semilla = -975;  				
