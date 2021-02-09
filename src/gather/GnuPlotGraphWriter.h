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
//! \file GnuPlotGraphWriter.h
//! \ingroup gather
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "DataProcessor.h"
#include "core/DataTypes.h"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>


namespace walberla {
namespace gather {


class GnuPlotGraphWriter : public DataProcessor
{
   public:
      typedef std::vector<real_t> DataPoint;

      GnuPlotGraphWriter(const std::string & filename, const std::string & fileEnding = "dat")
            : filename_(filename), fileEnding_(fileEnding), callNr(0) {}


      virtual void writeDataSet  (const std::vector<DataPoint> & dataset)
      {
         std::stringstream curFileName;
         curFileName << filename_ << '_' << callNr << "." << fileEnding_;

         std::ofstream outfile (curFileName.str().c_str() );
         for(unsigned int row = 0; row < dataset.size(); ++ row)
         {
            for(unsigned int col=0; col < dataset[row].size(); ++ col)
               outfile << dataset[row][col] << "\t";

            outfile << std::endl;
         }
         ++callNr;
      }

      void process(const std::vector<std::vector<real_t> > & data) override  { writeDataSet(data); }

   private:
      std::string filename_;
      std::string fileEnding_;
      uint_t      callNr;
};



} // namespace gather
} // namespace walberla


