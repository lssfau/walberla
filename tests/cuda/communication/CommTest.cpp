//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file
//! \author Paulo Carvalho <prcjunior@inf.ufpr.br>
//
//======================================================================================================================


#include "core/debug/TestSubsystem.h"
#include "core/Environment.h"
#include "core/mpi/Datatype.h"

#include "field/communication/MPIDatatypes.h"
#include "field/Field.h"

#include "cuda/GPUField.h"
#include "cuda/FieldCopy.h"

#define NUM_ITER  100
#define SIZE_X    16
#define SIZE_Y    16
#define SIZE_Z    16
#define LAYOUT    field::fzyx


using namespace walberla;

void hostToHost()
{
	Field<double, 1> hostField1(SIZE_X, SIZE_Y, SIZE_Z, 4321.0, LAYOUT);
	Field<double, 1> hostField2(SIZE_X, SIZE_Y, SIZE_Z, 0, LAYOUT);

	double startTime = MPI_Wtime();
	for (int i = 0; i < NUM_ITER; ++i) {
		hostField2.set(hostField1);
	}
	double endTime = MPI_Wtime();
	std::cout << __FUNCTION__ << ": " << endTime - startTime << std::endl;
}

void hostToDevice()
{
	Field<double, 1> hostField(SIZE_X, SIZE_Y, SIZE_Z, 4321.0, LAYOUT);
	cuda::GPUField<double> deviceField(SIZE_X, SIZE_Y, SIZE_Z, 1, 0, LAYOUT);

	double startTime = MPI_Wtime();
	for (int i = 0; i < NUM_ITER; ++i) {
		cuda::fieldCpy(deviceField, hostField);
	}
	double endTime = MPI_Wtime();
	std::cout << __FUNCTION__ << ": " << endTime - startTime << std::endl;
}

void deviceToHost()
{
	Field<double, 1> hostField(SIZE_X, SIZE_Y, SIZE_Z, 4321.0, LAYOUT);
	cuda::GPUField<double> deviceField(SIZE_X, SIZE_Y, SIZE_Z, 1, 0, LAYOUT);
	cuda::fieldCpy(deviceField, hostField);

	double startTime = MPI_Wtime();
	for (int i = 0; i < NUM_ITER; ++i) {
		cuda::fieldCpy(hostField, deviceField);
	}
	double endTime = MPI_Wtime();
	std::cout << __FUNCTION__ << ": " << endTime - startTime << std::endl;
}

void mpiHostToHost()
{
	Field<double, 1> hostField1(SIZE_X, SIZE_Y, SIZE_Z, 4321.0, LAYOUT);
	Field<double, 1> hostField2(SIZE_X, SIZE_Y, SIZE_Z, 0.0, LAYOUT);

	auto hostDatatype1 = mpi::Datatype ( field::communication::mpiDatatype( hostField1 ) );
	auto hostDatatype2 = mpi::Datatype ( field::communication::mpiDatatype( hostField2 ) );

	double startTime = MPI_Wtime();
	for (int i = 0; i < NUM_ITER; ++i) {
		MPI_Request request1;
		MPI_Isend( hostField1.data(), 1, hostDatatype1, 0, 0, MPI_COMM_WORLD, &request1 );

		MPI_Request request2;
		MPI_Irecv( hostField2.data(), 1, hostDatatype2, 0, 0, MPI_COMM_WORLD, &request2 );

		MPI_Wait( &request1, MPI_STATUS_IGNORE );
		MPI_Wait( &request2, MPI_STATUS_IGNORE );
	}
	double endTime = MPI_Wtime();
	std::cout << __FUNCTION__ << ": " << endTime - startTime << std::endl;
}

void mpiHostToDevice()
{
	Field<double, 1> hostField(SIZE_X, SIZE_Y, SIZE_Z, 4321.0, LAYOUT);
	cuda::GPUField<double> deviceField(SIZE_X, SIZE_Y, SIZE_Z, 1, 0, LAYOUT);

	auto hostDatatype = mpi::Datatype ( field::communication::mpiDatatype( hostField ) );
	auto deviceDatatype = mpi::Datatype ( field::communication::mpiDatatype( deviceField ) );

	double startTime = MPI_Wtime();
	for (int i = 0; i < NUM_ITER; ++i) {
		MPI_Request request1;
		MPI_Isend( hostField.data(), 1, hostDatatype, 0, 0, MPI_COMM_WORLD, &request1 );

		MPI_Request request2;
		MPI_Irecv( deviceField.data(), 1, deviceDatatype, 0, 0, MPI_COMM_WORLD, &request2 );

		MPI_Wait( &request1, MPI_STATUS_IGNORE );
		MPI_Wait( &request2, MPI_STATUS_IGNORE );
	}
	double endTime = MPI_Wtime();
	std::cout << __FUNCTION__ << ": " << endTime - startTime << std::endl;
}

void mpiDeviceToHost()
{
	Field<double, 1> hostField(SIZE_X, SIZE_Y, SIZE_Z, 4321.0, LAYOUT);
	cuda::GPUField<double> deviceField(SIZE_X, SIZE_Y, SIZE_Z, 1, 0, LAYOUT);

	auto hostDatatype = mpi::Datatype ( field::communication::mpiDatatype( hostField ) );
	auto deviceDatatype = mpi::Datatype ( field::communication::mpiDatatype( deviceField ) );

	double startTime = MPI_Wtime();
	for (int i = 0; i < NUM_ITER; ++i) {
		MPI_Request request1;
		MPI_Isend( deviceField.data(), 1, deviceDatatype, 0, 0, MPI_COMM_WORLD, &request1 );

		MPI_Request request2;
		MPI_Irecv( hostField.data(), 1, hostDatatype, 0, 0, MPI_COMM_WORLD, &request2 );

		MPI_Wait( &request1, MPI_STATUS_IGNORE );
		MPI_Wait( &request2, MPI_STATUS_IGNORE );
	}
	double endTime = MPI_Wtime();
	std::cout << __FUNCTION__ << ": " << endTime - startTime << std::endl;
}

void mpiDeviceToDevice()
{
	cuda::GPUField<double> deviceField1(SIZE_X, SIZE_Y, SIZE_Z, 1, 0, LAYOUT);
	cuda::GPUField<double> deviceField2(SIZE_X, SIZE_Y, SIZE_Z, 1, 0, LAYOUT);

	auto deviceDatatype1 = mpi::Datatype ( field::communication::mpiDatatype( deviceField1 ) );
	auto deviceDatatype2 = mpi::Datatype ( field::communication::mpiDatatype( deviceField2 ) );

	double startTime = MPI_Wtime();
	for (int i = 0; i < NUM_ITER; ++i) {
		MPI_Request request1;
		MPI_Isend( deviceField1.data(), 1, deviceDatatype1, 0, 0, MPI_COMM_WORLD, &request1 );

		MPI_Request request2;
		MPI_Irecv( deviceField2.data(), 1, deviceDatatype2, 0, 0, MPI_COMM_WORLD, &request2 );

		MPI_Wait( &request1, MPI_STATUS_IGNORE );
		MPI_Wait( &request2, MPI_STATUS_IGNORE );
	}
	double endTime = MPI_Wtime();
	std::cout << __FUNCTION__ << ": " << endTime - startTime << std::endl;
}

void mpiCopyHostToDevice()
{
	Field<double, 1> hostField1(SIZE_X, SIZE_Y, SIZE_Z, 4321.0, LAYOUT);
	Field<double, 1> hostField2(SIZE_X, SIZE_Y, SIZE_Z, 0.0, LAYOUT);
	cuda::GPUField<double> deviceField(SIZE_X, SIZE_Y, SIZE_Z, 1, 0, LAYOUT);

	auto hostDatatype1 = mpi::Datatype ( field::communication::mpiDatatype( hostField1 ) );
	auto hostDatatype2 = mpi::Datatype ( field::communication::mpiDatatype( hostField2 ) );

	double startTime = MPI_Wtime();
	for (int i = 0; i < NUM_ITER; ++i) {
		MPI_Request request1;
		MPI_Isend( hostField1.data(), 1, hostDatatype1, 0, 0, MPI_COMM_WORLD, &request1 );

		MPI_Request request2;
		MPI_Irecv( hostField2.data(), 1, hostDatatype2, 0, 0, MPI_COMM_WORLD, &request2 );

		MPI_Wait( &request1, MPI_STATUS_IGNORE );
		MPI_Wait( &request2, MPI_STATUS_IGNORE );

		cuda::fieldCpy(deviceField, hostField2);
	}
	double endTime = MPI_Wtime();
	std::cout << __FUNCTION__ << ": " << endTime - startTime << std::endl;
}

void mpiCopyDeviceToHost()
{
	Field<double, 1> hostField1(SIZE_X, SIZE_Y, SIZE_Z, 4321.0, LAYOUT);
	Field<double, 1> hostField2(SIZE_X, SIZE_Y, SIZE_Z, 0.0, LAYOUT);
	cuda::GPUField<double> deviceField(SIZE_X, SIZE_Y, SIZE_Z, 1, 0, LAYOUT);

	auto hostDatatype1 = mpi::Datatype ( field::communication::mpiDatatype( hostField1 ) );
	auto hostDatatype2 = mpi::Datatype ( field::communication::mpiDatatype( hostField2 ) );

	double startTime = MPI_Wtime();
	for (int i = 0; i < NUM_ITER; ++i) {

		MPI_Request request2;
		MPI_Irecv( hostField2.data(), 1, hostDatatype2, 0, 0, MPI_COMM_WORLD, &request2 );

		cuda::fieldCpy(hostField1, deviceField);

		MPI_Request request1;
		MPI_Isend( hostField1.data(), 1, hostDatatype1, 0, 0, MPI_COMM_WORLD, &request1 );

		MPI_Wait( &request1, MPI_STATUS_IGNORE );
		MPI_Wait( &request2, MPI_STATUS_IGNORE );
	}
	double endTime = MPI_Wtime();
	std::cout << __FUNCTION__ << ": " << endTime - startTime << std::endl;
}

int main( int argc, char ** argv )
{
   debug::enterTestMode();
   walberla::Environment walberlaEnv( argc, argv );

	WALBERLA_CHECK_EQUAL(MPIManager::instance()->numProcesses(), 2);

   hostToHost();
   hostToDevice();
   deviceToHost();
   mpiHostToHost();
   mpiHostToDevice();
   mpiDeviceToHost();
   mpiDeviceToDevice();
   mpiCopyHostToDevice();
   mpiCopyDeviceToHost();

   return 0;
}
