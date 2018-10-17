/*
 * @BEGIN LICENSE
 *
 * embedding by Psi4 Developer, a plugin to:
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

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.hpp"
#include "emb.h"

namespace psi{ namespace embedding {

extern "C" PSI_API
int read_options(std::string name, Options& options)
{
    if (name == "EMBEDDING"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
		options.add_str("CUTOFF_BY", "THRESHOLD");
		options.add_int("NUM_OCC", 0);
		options.add_int("NUM_VIR", 0);
		options.add_double("THRESHOLD", 0.5);
		options.add_str("REFERENCE", "HF");
		options.add_bool("WRITE_FREEZE_MO", true);
    }

    return true;
}

extern "C" PSI_API
SharedWavefunction embedding(SharedWavefunction ref_wfn, Options& options)
{
	emb e1(ref_wfn, options);
	e1.compute_energy();
	return ref_wfn;
}

}} // End namespaces

