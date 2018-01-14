/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2018 Sebastian Wouters

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef DIIS_CHEMPS2_H
#define DIIS_CHEMPS2_H

#include "Options.h"

namespace CheMPS2{
/** DIIS class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date July 8, 2014
    
    The DIIS class calculates the new step based on a direct inversion of iterative subspaces (DIIS). Info can be found at
    - http://vergil.chemistry.gatech.edu/notes/diis/node3.html
    - Yanai, Kurashige, Ghosh, Chan, International Journal of Quantum Chemistry 109, 2178-2190 (2009) [ http://dx.doi.org/10.1002/qua.22099 ].
*/
   class DIIS{

      public:
      
         //! Constructor
         /** \param numVarsParamIn The number of variables in the parameter vectors
             \param numVarsErrorIn The number of variables in the error vectors
             \param numVecsIn The max. number of parameter and error vectors */
         DIIS(const int numVarsParamIn, const int numVarsErrorIn, const int numVecsIn);
         
         //! Destructor
         virtual ~DIIS();
         
         //! Get the number of variables in the parameter vectors
         /** \return The number of variables in the parameter vectors */
         int getNumVarsParam() const;
         
         //! Get the number of variables in the error vectors
         /** \return The number of variables in the error vectors */
         int getNumVarsError() const;
         
         //! Get the max. number of parameter and error vectors which are kept
         /** \return The max. number of parameter and error vectors which are kept */
         int getNumVecs() const;
         
         //! Get the current number of vectors which are used for DIIS
         /** \return The current number of vectors which are used for DIIS */
         int getCurrentNumVecs() const;
         
         //! Get pointer to the last linco of parameter vectors
         /** \return Pointer to the last linco of parameter vectors */
         double * getLastLinco();
         
         //! Append a new error and parameter vector
         /** \param newError The new error vector which should be appended
             \param newParam The new parameter vector which should be appended */
         void appendNew(double * newError, double * newParam);
         
         //! Calculate the new parameter vector, based on the just appended error and parameter vectors.
         /** \param newParam The new parameter vector */
         void calculateParam(double * newParam);
         
         //! Save the DIIS object to disk
         /** \param filename Filename to store the DIIS checkpoint to */
         void saveDIIS(const string filename=DMRGSCF_diis_storage_name) const;
         
         //! Load the DIIS object from disk
         /** \param filename Filename to load the DIIS checkpoint from */
         void loadDIIS(const string filename=DMRGSCF_diis_storage_name);

      private:

         //The number of variables in the parameter vectors
         int numVarsParam;
         
         //The number of variables in the error vectors
         int numVarsError;
         
         //The max. number of vectors which should be kept
         int numVecs;
         
         //The current number of vectors in this object
         int currentNumVecs;
      
         //The error vectors
         double ** errorVectors;
         
         //The parameter vectors
         double ** paramVectors;
         
         //The last linco of parameter vectors
         double * lastLinco;
         
   };
}

#endif
