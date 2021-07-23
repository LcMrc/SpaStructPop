#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Initialization(int WildTypeVector[], int MutantVector[], int InitialNumberOfWildTypes, int InitialNumberOfMutants, double *Time)
{
	WildTypeVector[0] = InitialNumberOfWildTypes;
	WildTypeVector[1] = InitialNumberOfWildTypes;
	WildTypeVector[2] = InitialNumberOfWildTypes;
	WildTypeVector[3] = InitialNumberOfWildTypes;
	WildTypeVector[4] = 0;
	MutantVector[0] = 0;
	MutantVector[1] = 0;
	MutantVector[2] = 0;
	MutantVector[3] = 0;
	MutantVector[4] = InitialNumberOfMutants;	
	*Time = 0;
}

int Sum(int PopulationVector[], int NumberOfDemes, int *TotalNumberOfIndividuals)
{

	int i, q = 0;

	for (i = 0; i < NumberOfDemes; i++)
	{ 
		q = q+PopulationVector[i];
	}

	*TotalNumberOfIndividuals = q;

}

double ComputeTotalDivisionRate(double Fitness, int WildTypeVector[], int MutantVector[], int CarryingCapacity, int NumberOfDemes)
{
	double TotalDivisionRate = 0;
	int i;

	for (i = 0; i < NumberOfDemes; i++)
	{
		TotalDivisionRate = TotalDivisionRate+Fitness*(1-(double)(WildTypeVector[i]+MutantVector[i])/CarryingCapacity)*WildTypeVector[i];
	}

	return TotalDivisionRate;
}

double ComputeTotalDeathRate(double DeathRate, int PopulationVector[], int NumberOfDemes)
{
	double TotalDeathRate = 0;
	int i;

	for (i = 0; i < NumberOfDemes; i++)
	{
		TotalDeathRate =TotalDeathRate+DeathRate*PopulationVector[i];
	}

	return TotalDeathRate;
}

double ComputeTotalMigrationRate(int WildTypeVector[], int MutantVector[], double MigrationRate, int NumberOfDemes)
{
	double TotalMigrationRate = 0;

	TotalMigrationRate = MigrationRate*(NumberOfDemes-1)*(WildTypeVector[0]+MutantVector[0]+WildTypeVector[1]+MutantVector[1]+WildTypeVector[2]+MutantVector[2]+WildTypeVector[3]+MutantVector[3]+WildTypeVector[4]+MutantVector[4]);

	return TotalMigrationRate;
}

void ComputeTotalTransitionRate(double WildTypeFitness, double WildTypeDeathRate, int WildTypeVector[], double MutantFitness, double MutantDeathRate, int MutantVector[], int CarryingCapacity, int NumberOfDemes, double MigrationRate, double *TransitionRatesVector)
{
	TransitionRatesVector[0] = ComputeTotalDivisionRate(WildTypeFitness, WildTypeVector, MutantVector, CarryingCapacity, NumberOfDemes);
	TransitionRatesVector[1] = ComputeTotalDeathRate(WildTypeDeathRate, WildTypeVector, NumberOfDemes);
	TransitionRatesVector[2] = ComputeTotalDivisionRate(MutantFitness, MutantVector, WildTypeVector, CarryingCapacity, NumberOfDemes);
	TransitionRatesVector[3] = ComputeTotalDeathRate(MutantDeathRate, MutantVector, NumberOfDemes);
	TransitionRatesVector[4] = ComputeTotalMigrationRate(WildTypeVector, MutantVector, MigrationRate, NumberOfDemes);
}

void ComputeCumulativeSum(double *TransitionRatesVector, int VectorSize) 
{
    if(VectorSize <= 0) return;
    ComputeCumulativeSum(TransitionRatesVector, VectorSize-1);
    TransitionRatesVector[VectorSize] += TransitionRatesVector[VectorSize-1];
}

double ComputeRandomNumber()
{
    return (double)rand()/(double)RAND_MAX;
}

int SamplingEvent(double *TransitionRatesVector, int VectorSize) 
{
	int i = 0;
	double RandomNumber1, RandomNumber2;
	
	RandomNumber1 = ComputeRandomNumber();
	RandomNumber2 = RandomNumber1*TransitionRatesVector[VectorSize-1];

	while (TransitionRatesVector[i] < RandomNumber2)
	{
		i++;
	}

	return i;
}

double UpdateTime(double *Time, double TotalTransitionRate)
{
	double RandomNumber, PopulationVector;
	RandomNumber = ComputeRandomNumber();
	PopulationVector = 1/RandomNumber;
	*Time += 1/TotalTransitionRate*log(PopulationVector);
}

void PartialTransitionRateBis(double *TransitionRatesVector, int *WildTypeVector, int *MutantVector, int CarryingCapacity, int NumberOfDemes, int Event)
{
	int i;

	if (Event == 0)
	{
		for (i = 0; i < NumberOfDemes; i++)
		{
			TransitionRatesVector[i] = (1-(double)(WildTypeVector[i]+MutantVector[i])/CarryingCapacity)*WildTypeVector[i];
		}
	}
	else if (Event == 1)
	{
		for (i = 0; i < NumberOfDemes; i++)
		{
			TransitionRatesVector[i] = WildTypeVector[i];
		}
	}
	else if (Event == 2)
	{
		for (i = 0; i < NumberOfDemes; i++)
		{
			TransitionRatesVector[i] = (1-(double)(WildTypeVector[i]+MutantVector[i])/CarryingCapacity)*MutantVector[i];
		}
	}
	else
	{
		for (i = 0; i < NumberOfDemes; i++)
		{
			TransitionRatesVector[i] = MutantVector[i];
		}
	}
}

void ComputeMigrationRates(double *TransitionRatesVector, int *WildTypeVector, int *MutantVector, double MigrationRate)
{
	TransitionRatesVector[0] = MigrationRate*(WildTypeVector[0]+MutantVector[0]);
	TransitionRatesVector[1] = MigrationRate*(WildTypeVector[0]+MutantVector[0]);
	TransitionRatesVector[2] = MigrationRate*(WildTypeVector[0]+MutantVector[0]);
	TransitionRatesVector[3] = MigrationRate*(WildTypeVector[0]+MutantVector[0]);
	TransitionRatesVector[4] = MigrationRate*(WildTypeVector[1]+MutantVector[1]);
	TransitionRatesVector[5] = MigrationRate*(WildTypeVector[1]+MutantVector[1]);
	TransitionRatesVector[6] = MigrationRate*(WildTypeVector[1]+MutantVector[1]);
	TransitionRatesVector[7] = MigrationRate*(WildTypeVector[1]+MutantVector[1]);	
	TransitionRatesVector[8] = MigrationRate*(WildTypeVector[2]+MutantVector[2]);
	TransitionRatesVector[9] = MigrationRate*(WildTypeVector[2]+MutantVector[2]);
	TransitionRatesVector[10] = MigrationRate*(WildTypeVector[2]+MutantVector[2]);
	TransitionRatesVector[11] = MigrationRate*(WildTypeVector[2]+MutantVector[2]);
	TransitionRatesVector[12] = MigrationRate*(WildTypeVector[3]+MutantVector[3]);
	TransitionRatesVector[13] = MigrationRate*(WildTypeVector[3]+MutantVector[3]);
	TransitionRatesVector[14] = MigrationRate*(WildTypeVector[3]+MutantVector[3]);
	TransitionRatesVector[15] = MigrationRate*(WildTypeVector[3]+MutantVector[3]);
	TransitionRatesVector[16] = MigrationRate*(WildTypeVector[4]+MutantVector[4]);
	TransitionRatesVector[17] = MigrationRate*(WildTypeVector[4]+MutantVector[4]);
	TransitionRatesVector[18] = MigrationRate*(WildTypeVector[4]+MutantVector[4]);
	TransitionRatesVector[19] = MigrationRate*(WildTypeVector[4]+MutantVector[4]);
}

void ExecuteDivisionOrDeath(int *WildTypeVector, int *MutantVector, int Event, int TargetedDeme)
{
	if (Event == 0)
	{
		WildTypeVector[TargetedDeme] += 1;
	}
	else if (Event == 1)
	{
		WildTypeVector[TargetedDeme] -= 1;
	}
	else if (Event == 2)
	{
		MutantVector[TargetedDeme] += 1;
	}
	else
	{
		MutantVector[TargetedDeme] -= 1;
	}
}

void ExecuteMigration(int *WildTypeVector, int *MutantVector, int Event)
{
	int DepartureDeme, ArrivalDeme;
	double RandomNumber;

	if (Event == 0)
	{
		DepartureDeme = 0;
		ArrivalDeme = 1;
	}
	else if (Event == 1)
	{
		DepartureDeme = 0;
		ArrivalDeme = 2;
	}
	else if (Event == 2)
	{
		DepartureDeme = 0;
		ArrivalDeme = 3;
	}
	else if (Event == 3)
	{
		DepartureDeme = 0;
		ArrivalDeme = 4;
	}
	else if (Event == 4)
	{
		DepartureDeme = 1;
		ArrivalDeme = 0;
	}
	else if (Event == 5)
	{
		DepartureDeme = 1;
		ArrivalDeme = 2;
	}
	else if (Event == 6)
	{
		DepartureDeme = 1;
		ArrivalDeme = 3;
	}
	else if (Event == 7)
	{
		DepartureDeme = 1;
		ArrivalDeme = 4;
	}
	else if (Event == 8)
	{
		DepartureDeme = 2;
		ArrivalDeme = 0;
	}
	else if (Event == 9)
	{
		DepartureDeme = 2;
		ArrivalDeme = 1;
	}
	else if (Event == 10)
	{
		DepartureDeme = 2;
		ArrivalDeme = 3;
	}
	else if (Event == 11)
	{
		DepartureDeme = 2;
		ArrivalDeme = 4;
	}
	else if (Event == 12)
	{
		DepartureDeme = 3;
		ArrivalDeme = 0;
	}
	else if (Event == 13)
	{
		DepartureDeme = 3;
		ArrivalDeme = 1;
	}
	else if (Event == 14)
	{
		DepartureDeme = 3;
		ArrivalDeme = 2;
	}
	else if (Event == 15)
	{
		DepartureDeme = 3;
		ArrivalDeme = 4;
	}
	else if (Event == 16)
	{
		DepartureDeme = 4;
		ArrivalDeme = 0;
	}
	else if (Event == 17)
	{
		DepartureDeme = 4;
		ArrivalDeme = 1;
	}
	else if (Event == 18)
	{
		DepartureDeme = 4;
		ArrivalDeme = 2;
	}
	else if (Event == 19)
	{
		DepartureDeme = 4;
		ArrivalDeme = 3;
	}

	RandomNumber = ComputeRandomNumber();

	if (RandomNumber <= WildTypeVector[DepartureDeme]/(WildTypeVector[DepartureDeme]+MutantVector[DepartureDeme]))
	{
		WildTypeVector[DepartureDeme] -= 1;
		WildTypeVector[ArrivalDeme] += 1;
	}
	else
	{
		MutantVector[DepartureDeme] -= 1;
		MutantVector[ArrivalDeme] += 1;
	}
}

int main()
{

	const double WildTypeFitness = 1.0, WildTypeDeathRate = 0.1, MutantFitness = 1.01, MutantDeathRate = 0.1, MigrationRate = 0.000001;
	const int CarryingCapacity = 100, InitialNumberOfWildTypes = 90, InitialNumberOfMutants = 90, NumberOfDemes = 5, NumberOfRealizations = 1000;
	double Time = 0, TotalTransitionRateVector[5], DivisionAndDeathTransitionRateVector[5], MigrationRateVector[20], FinalTimeVector[1000];
	int WildTypeVector[5] = {InitialNumberOfWildTypes, InitialNumberOfWildTypes, InitialNumberOfWildTypes, InitialNumberOfWildTypes, 0}, MutantVector[5] = {0, 0, 0, 0, InitialNumberOfMutants}, k, i, j, Event, TargetedDeme, TotalNumberOfWildTypes = 0, TotalNumberOfMutants = 0, FinalMutantNumberVector[1000];

	for (k = 0; k < NumberOfRealizations; k++)
	{
		Initialization(WildTypeVector, MutantVector, InitialNumberOfWildTypes, InitialNumberOfMutants, &Time);
		Sum(WildTypeVector, NumberOfDemes, &TotalNumberOfWildTypes);
		Sum(MutantVector, NumberOfDemes, &TotalNumberOfMutants);

		while (TotalNumberOfWildTypes != 0 && TotalNumberOfMutants != 0)
		{
			ComputeTotalTransitionRate(WildTypeFitness, WildTypeDeathRate, WildTypeVector, MutantFitness, MutantDeathRate, MutantVector, CarryingCapacity, NumberOfDemes, MigrationRate, TotalTransitionRateVector);
			ComputeCumulativeSum(TotalTransitionRateVector, NumberOfDemes-1);
			Event = SamplingEvent(TotalTransitionRateVector, NumberOfDemes); 

			if (Event < 4)
			{
				PartialTransitionRateBis(DivisionAndDeathTransitionRateVector, WildTypeVector, MutantVector, CarryingCapacity, NumberOfDemes, Event);
				ComputeCumulativeSum(DivisionAndDeathTransitionRateVector, NumberOfDemes-1);
				TargetedDeme = SamplingEvent(DivisionAndDeathTransitionRateVector, NumberOfDemes);
				ExecuteDivisionOrDeath(WildTypeVector, MutantVector, Event, TargetedDeme);
			}
			else
			{
				ComputeMigrationRates(MigrationRateVector, WildTypeVector, MutantVector, MigrationRate);
				ComputeCumulativeSum(MigrationRateVector, 20-1);
				TargetedDeme = SamplingEvent(MigrationRateVector, 20);
				ExecuteMigration(WildTypeVector, MutantVector, TargetedDeme);
			}

			UpdateTime(&Time, TotalTransitionRateVector[NumberOfDemes-1]);

			Sum(WildTypeVector, NumberOfDemes, &TotalNumberOfWildTypes);
			Sum(MutantVector, NumberOfDemes, &TotalNumberOfMutants);	 
		}
		
		FinalMutantNumberVector[k] = MutantVector[4];
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
