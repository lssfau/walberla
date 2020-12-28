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
//! \file TimeloopAndSweepRegister.cpp
//! \ingroup timeloop
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief test cases that test the registering of Sweeps at timeloop
//
//======================================================================================================================

#include "timeloop/SweepTimeloop.h"
#include "core/DataTypes.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"

#include <string>
#include <vector>


using namespace std;
using namespace walberla;


class GeneralSweep
{
   public:
      GeneralSweep(const string & name, vector<string> & vec)
         : myName(name), outVec(vec)
      {}


      void operator() (IBlock * ){
         outVec.push_back(myName);
      }

   private:
      string myName;
      vector<string> & outVec;
};


class GeneralFunction
{
   public:
      GeneralFunction(const string & name, vector<string> & vec)
      : myName(name), outVec(vec)
      {}

      void operator() (){
         outVec.push_back(myName);
      }

   private:
      string myName;
      vector<string> & outVec;
};



int main()
{
   debug::enterTestMode();

   vector<string> expectedSequence;
   vector<string> sequence;

   SUID cpuSelect("CPU");
   SUID gpuSelect("GPU");

   SUID srtSelect("SRT");
   SUID mrtSelect("MRT");


   //FIXME put a real test in here
   shared_ptr<SweepTimeloop> tl = make_shared<SweepTimeloop>(shared_ptr<BlockStorage>(),100);

   typedef SweepTimeloop::SweepHandle SH;

   SH sweep1 =  tl->addSweep(        GeneralSweep("CPU1",sequence), cpuSelect);
                tl->addSweep(sweep1, GeneralSweep("GPU1",sequence), gpuSelect);

   SH sweep2 =  tl->addSweep(        GeneralSweep("CPU2",sequence), cpuSelect);
                tl->addSweep(sweep2, GeneralSweep("GPU2",sequence), gpuSelect);


   tl->addFuncBeforeSweep(sweep1,GeneralFunction("Pre1",sequence));
   tl->addFuncAfterSweep (sweep1,GeneralFunction("Post1",sequence));

   tl->addFuncBeforeSweep(sweep2,GeneralFunction("Pre2",sequence));
   tl->addFuncAfterSweep (sweep2,GeneralFunction("Post2",sequence));

   typedef Timeloop::FctHandle FH;
   FH preTs  = tl->addFuncBeforeTimeStep(       GeneralFunction("PreTimestepCPU",sequence),cpuSelect,srtSelect);
               tl->addFuncBeforeTimeStep(preTs, GeneralFunction("PreTimestepGPU",sequence),gpuSelect,srtSelect);
   FH postTs = tl->addFuncAfterTimeStep (       GeneralFunction("PostTimestepCPU",sequence),cpuSelect,mrtSelect);
               tl->addFuncAfterTimeStep (postTs,GeneralFunction("PostTimestepGPU",sequence),gpuSelect,mrtSelect);


   //----------  First Run - CPU Selector ---------------------------------------------

   expectedSequence.push_back("PreTimestepCPU");
   expectedSequence.push_back("Pre1");
   expectedSequence.push_back("CPU1");
   expectedSequence.push_back("Post1");
   expectedSequence.push_back("Pre2");
   expectedSequence.push_back("CPU2");
   expectedSequence.push_back("Post2");
   expectedSequence.push_back("PostTimestepCPU");

   tl->singleStep(cpuSelect);

   WALBERLA_CHECK_EQUAL(expectedSequence.size(), sequence.size());
   WALBERLA_CHECK( equal(expectedSequence.begin(),expectedSequence.end(), sequence.begin() ) );

   expectedSequence.clear();
   sequence.clear();


   // ------------ Second Run - GPU Selector -------------------------------------------

   expectedSequence.push_back("PreTimestepGPU");
   expectedSequence.push_back("Pre1");
   expectedSequence.push_back("GPU1");
   expectedSequence.push_back("Post1");
   expectedSequence.push_back("Pre2");
   expectedSequence.push_back("GPU2");
   expectedSequence.push_back("Post2");
   expectedSequence.push_back("PostTimestepGPU");

   tl->singleStep(gpuSelect);

   WALBERLA_CHECK_EQUAL(expectedSequence.size(), sequence.size());
   WALBERLA_CHECK( equal(expectedSequence.begin(),expectedSequence.end(), sequence.begin() ) );

   expectedSequence.clear();
   sequence.clear();


   // ------------ Second Run - GPU and SRT    -------------------------------------------


   expectedSequence.push_back("Pre1");
   expectedSequence.push_back("GPU1");
   expectedSequence.push_back("Post1");
   expectedSequence.push_back("Pre2");
   expectedSequence.push_back("GPU2");
   expectedSequence.push_back("Post2");
   expectedSequence.push_back("PostTimestepGPU");

   tl->singleStep(gpuSelect + srtSelect);

   WALBERLA_CHECK_EQUAL(expectedSequence.size(), sequence.size());
   WALBERLA_CHECK( equal(expectedSequence.begin(),expectedSequence.end(), sequence.begin() ) );

   expectedSequence.clear();
   sequence.clear();


   return 0;
}



