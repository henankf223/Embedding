/*
* @BEGIN LICENSE
*
* new_ivocalc by Psi4 Developer, a plugin to:
*
* Psi4: an open-source quantum chemistry software package
*
* Copyright (c) 2007-2018 The Psi4 Developers.
*
* The copyrights for code used from other parties are included in
* the corresponding files.
*
* This file is part of Psi4.
*
* Psi4 is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, version 3.
*
* Psi4 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License along
* with Psi4; if not, write to the Free Software Foundation, Inc.,
* 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*
* @END LICENSE
*/

#include "psi4/libmints/wavefunction.h"

namespace psi {
	namespace embedding {

		/**
		* @The IVO class
		* This class implements Improved Virtual Orbitals with one hole
		*/
		class emb : public Wavefunction {
		public:
			// ==> Class Constructor and Destructor <==

			/**
			* Constructor
			* @param ref_wfn The reference wavefunction object
			* @param options The main options object
			*/
			emb(SharedWavefunction ref_wfn, Options& options);

			/// Destructor
			~emb();

			// ==> Class Interface <==

			/// Run calculation
			double compute_energy();

		private:
			// ==> Class data <==
			int print;
			int num_occ;
			int num_vir;
			double thresh;

			// ==> Class functions <==

		};

	}
}
