#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Initialization(int WildTypeVector[], int MutantVector[], int InitialNumberOfWildTypes, int InitialNumberOfMutants, double *Time)
{
	WildTypeVector[0] = InitialNumberOfWildTypes;
	WildTypeVector[1] = 0;
	MutantVector[0] = 0;
	MutantVector[1] = InitialNumberOfMutants;	
	*Time = 0;
}

int ComputeSum(int PopulationVector[], int NumberOfDemes, int *Sum)
{

	int i, q = 0;

	for (i = 0; i < NumberOfDemes; i++)
	{ 
		q = q+PopulationVector[i];
	}

	*Sum = q;

}

void TotalTransitionRate(double WildTypeFitness, double WildTypeDeathRate, int WildTypeVector[], double MutantFitness, double MutantDeathRate, int MutantVector[], int LargeCarryingCapacity, int SmallCarryingCapacity, int NumberOfDemes, double SmallToLargeMigrationRate, double LargeToSmallMigrationRate, double *TransitionRateVector)
{
	TransitionRateVector[0] = WildTypeFitness*(1-(double)(WildTypeVector[0]+MutantVector[0])/LargeCarryingCapacity)*WildTypeVector[0];
	TransitionRateVector[1] = WildTypeDeathRate*WildTypeVector[0];
	TransitionRateVector[2] = WildTypeFitness*(1-(double)(WildTypeVector[1]+MutantVector[1])/SmallCarryingCapacity)*WildTypeVector[1];
	TransitionRateVector[3] = WildTypeDeathRate*WildTypeVector[1];
	TransitionRateVector[4] = MutantFitness*(1-(double)(WildTypeVector[0]+MutantVector[0])/LargeCarryingCapacity)*MutantVector[0];
	TransitionRateVector[5] = MutantDeathRate*MutantVector[0];
	TransitionRateVector[6] = MutantFitness*(1-(double)(WildTypeVector[1]+MutantVector[1])/SmallCarryingCapacity)*MutantVector[1];
	TransitionRateVector[7] = MutantDeathRate*MutantVector[1];
	TransitionRateVector[8] = LargeToSmallMigrationRate*WildTypeVector[0];
	TransitionRateVector[9] = SmallToLargeMigrationRate*WildTypeVector[1];
	TransitionRateVector[10] = LargeToSmallMigrationRate*MutantVector[0];
	TransitionRateVector[11] = SmallToLargeMigrationRate*MutantVector[1];
}

void ComputeCumulativeSum(double *TransitionRateVector, int NumberOfEvents) 
{
    if(NumberOfEvents <= 0) return;
    ComputeCumulativeSum(TransitionRateVector, NumberOfEvents-1);
    TransitionRateVector[NumberOfEvents] += TransitionRateVector[NumberOfEvents-1];
}

double ComputeRandomNumber()
{
    return (double)rand()/(double)RAND_MAX;
}

int SamplingEvent(double *TransitionRateVector, int NumberOfEvents) 
{
	int i = 0;
	double RandomNumber1, RandomNumber2;
	
	RandomNumber1 = ComputeRandomNumber();
	RandomNumber2 = RandomNumber1*TransitionRateVector[NumberOfEvents-1];

	while (TransitionRateVector[i] < RandomNumber2)
	{
		i++;
	}

	return i;
}

double UpdateTime(double *Time, double TotalTransitionRate)
{
	double RandomNumber;
	RandomNumber = ComputeRandomNumber();
	*Time += 1/TotalTransitionRate*log(1/RandomNumber);
}

void ExecuteReaction(int *WildTypeVector, int *MutantVector, int Reaction)
{
	if (Reaction == 0)
	{
		WildTypeVector[0] += 1;
	}
	else if (Reaction == 1)
	{
		WildTypeVector[0] -= 1;
	}
	else if (Reaction == 2)
	{
		WildTypeVector[1] += 1;
	}
	else if (Reaction == 3)
	{
		WildTypeVector[1] -= 1;
	}
	else if (Reaction == 4)
	{
		MutantVector[0] += 1;
	}
	else if (Reaction == 5)
	{
		MutantVector[0] -= 1;
	}
	else if (Reaction == 6)
	{
		MutantVector[1] += 1;
	}
	else if (Reaction == 7)
	{
		MutantVector[1] -= 1;
	}
	else if (Reaction == 8)
	{
		WildTypeVector[0] -= 1;
		WildTypeVector[1] += 1;	
	}
	else if (Reaction == 9)
	{
		WildTypeVector[0] += 1;
		WildTypeVector[1] -= 1;
	}
	else if (Reaction == 10)
	{
		MutantVector[0] -= 1;
		MutantVector[1] += 1;
	}
	else if (Reaction == 11)
	{
		MutantVector[0] += 1;
		MutantVector[1] -= 1;
	}
}

int main()
{

	const double WildTypeFitness = 1.0, WildTypeDeathRate = 0.1, MutantFitness = 1.01, MutantDeathRate = 0.1, SmallToLargeMigrationRate = 0.000008, LargeToSmallMigrationRate = 0.000001;
	const int LargeCarryingCapacity = 400, SmallCarryingCapacity = 100, InitialNumberOfWildTypes = 360, InitialNumberOfMutants = 90, NumberOfDemes = 2, NumberOfRealizations = 1000;
	double Time = 0, TransitionRateVector[12], FinalTimeVector[1000];
	int WildTypeVector[2] = {InitialNumberOfWildTypes, 0}, MutantVector[2] = {0, InitialNumberOfMutants}, k, i, j, Reaction, ir2, TotalNumberOfWildTypes = 0, TotalNumberOfMutants = 0, FinalMutantNumberVector[1000];

	for (k = 0; k < NumberOfRealizations; k++)
	{
		Initialization(WildTypeVector, MutantVector, InitialNumberOfWildTypes, InitialNumberOfMutants, &Time);
		ComputeSum(WildTypeVector, NumberOfDemes, &TotalNumberOfWildTypes);
		ComputeSum(MutantVector, NumberOfDemes, &TotalNumberOfMutants);

		while (TotalNumberOfWildTypes != 0 && TotalNumberOfMutants != 0)
		{
			TotalTransitionRate(WildTypeFitness, WildTypeDeathRate, WildTypeVector, MutantFitness, MutantDeathRate, MutantVector, LargeCarryingCapacity, SmallCarryingCapacity, NumberOfDemes, SmallToLargeMigrationRate, LargeToSmallMigrationRate, TransitionRateVector);
			ComputeCumulativeSum(TransitionRateVector, 11);
			Reaction = SamplingEvent(TransitionRateVector, 12); 
			ExecuteReaction(WildTypeVector, MutantVector, Reaction);

			UpdateTime(&Time, TransitionRateVector[11]);

			ComputeSum(WildTypeVector, NumberOfDemes, &TotalNumberOfWildTypes);
			ComputeSum(MutantVector, NumberOfDemes, &TotalNumberOfMutants);	 
		}
		
		FinalMutantNumberVector[k] = MutantVector[1];
		FinalTimeVector[k] = Time;
	}

	FILE *fp1;
	fp1 = fopen("FinalMutantNumberVector.txt", "w");
	for (i = 0; i < NumberOfRealizations; i++)
	{
		fprintf(fp1, "%d\n", FinalMutantNumberVector[i]);
	}
	fclose(fp1); 

	FILE *fp2;
	fp2 = fopen("FinalTimeVector.txt", "w");
	for (i = 0; i < NumberOfRealizations; i++)
	{
		fprintf(fp2, "%f\n", FinalTimeVector[i]);
	}
	fclose(fp2);

	return 0;
} 
