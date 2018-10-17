/*
 * @BEGIN LICENSE
 *
 * ivocalc by Psi4 Developer, a plugin to:
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
#include "psi4/libmints/vector.h" 
#include "psi4/libmints/matrix.h" 
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h" 
#include "psi4/libmints/wavefunction.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "emb.h"

namespace psi{ namespace embedding {

emb::emb(SharedWavefunction ref_wfn, Options& options): Wavefunction(options) {
	outfile->Printf("\n ------ Orbital Localization and Embedding ------ \n");
	shallow_copy(ref_wfn);
	reference_wavefunction_ = ref_wfn;

	//1. Get necessary information
	print = options_.get_int("PRINT");
	thresh = options_.get_double("THRESHOLD");
	num_occ = options_.get_int("NUM_OCC");
	num_vir = options_.get_int("NUM_VIR");

	std::shared_ptr<PSIO> psio(_default_psio_lib_);
	if (!ref_wfn)
		throw PSIEXCEPTION("SCF has not been run yet!");
}

emb::~emb() {}

double emb::compute_energy() {
	std::shared_ptr<Molecule> mol = reference_wavefunction_->molecule();
	int nfrag = mol->nfragments();

	if (nfrag == 1) {
		outfile->Printf("Warning! A input molecule with fragments (--- in atom list) is required "
			"for embedding!");
	}

	outfile->Printf(
		"\n The input molecule have %d fragments, assigning the first fragment as system! \n",
		nfrag);

	std::vector<int> none_list = {};
	std::vector<int> sys_list = { 0 };
	std::vector<int> env_list = { 1 }; // change to 1-nfrag in the future!

	std::shared_ptr<Molecule> mol_sys = mol->extract_subsets(sys_list, none_list);
	std::shared_ptr<Molecule> mol_env = mol->extract_subsets(env_list, none_list);
	outfile->Printf("\n System Fragment \n");
	mol_sys->print();
	outfile->Printf("\n Environment Fragment(s) \n");
	mol_env->print();

	std::shared_ptr<BasisSet> basis = reference_wavefunction_->basisset();
	Dimension nmopi = reference_wavefunction_->nmopi();
	int nirrep = reference_wavefunction_->nirrep();
	int nbf = basis->nbf();
	outfile->Printf("\n number of basis on all atoms: %d", nbf);

	int natom_sys = mol_sys->natom();
	int count_basis = 0;
	for (int mu = 0; mu < nbf; mu++) {
		int A = basisset_->function_to_center(mu);
		outfile->Printf("\n  Function %d is on atom %d", mu, A);
		if (A < natom_sys) {
			count_basis += 1;
		}
	}
	outfile->Printf("\n number of basis on \"system\" atoms: %d", count_basis);

	//Naive partition
	SharedMatrix S_ao = reference_wavefunction_->S();
	Dimension noccpi = reference_wavefunction_->doccpi();
	Dimension zeropi = nmopi - nmopi;
	Dimension nvirpi = nmopi - noccpi;
	Dimension sys_mo = nmopi;
	sys_mo[0] = count_basis;
	Slice sys(zeropi, sys_mo);
	Slice env(sys_mo, nmopi);
	Slice allmo(zeropi, nmopi);
	Slice occ(zeropi, noccpi);
	Slice vir(noccpi, nmopi);
	SharedMatrix S_sys = S_ao->get_block(sys, sys);
	SharedMatrix L(new Matrix("L", nirrep, sys_mo, sys_mo));
	SharedVector lm(new Vector("lambda", nirrep, sys_mo));
	SharedVector lminvhalf(new Vector("lambda inv half", nirrep, sys_mo));
	SharedMatrix LM(new Matrix("LM", nirrep, sys_mo, sys_mo));

	//Construct S_sys^-1/2
	S_sys->diagonalize(L, lm);
	for (int i = 0; i < sys_mo[0]; ++i) {
		double tmp = 1.0 / lm->get(0, i);
		lminvhalf->set(0, i, tmp);
	}
	LM->set_diagonal(lminvhalf);
	SharedMatrix S_sys_invhalf = Matrix::triplet(L, LM, L, false, false, true);

	SharedMatrix S_sys_in_all(new Matrix("S system in fullsize", nirrep, nmopi, nmopi));
	S_sys_in_all->set_block(sys, sys, S_sys_invhalf);
	//outfile->Printf("\n S_sys in large \n");
	//S_sys_in_all->print();

	//Build P_pq
	SharedMatrix Ca_t = reference_wavefunction_->Ca();
	S_sys_in_all->transform(S_ao);
	S_sys_in_all->transform(Ca_t);

	//Diagonalize P_pq for occ and vir part, respectively.
	SharedMatrix P_oo = S_sys_in_all->get_block(occ, occ);
	SharedMatrix Uo(new Matrix("Uo", nirrep, noccpi, noccpi));
	SharedVector lo(new Vector("lo", nirrep, noccpi));
	P_oo->diagonalize(Uo, lo, descending);
	lo->print();

	SharedMatrix P_vv = S_sys_in_all->get_block(vir, vir);
	SharedMatrix Uv(new Matrix("Uv", nirrep, nvirpi, nvirpi));
	SharedVector lv(new Vector("lv", nirrep, nvirpi));
	P_vv->diagonalize(Uv, lv, descending);
	lv->print();

	SharedMatrix U_all(new Matrix("U with Pab", nirrep, nmopi, nmopi));
	U_all->set_block(occ, occ, Uo);
	U_all->set_block(vir, vir, Uv);

	//Update MO coeffs
	Ca_->copy(Matrix::doublet(Ca_t, U_all, false, false));
	//Based on threshold or num_occ/num_vir, decide the mo_space_info
	std::vector<int> index_trace_occ = {};
	std::vector<int> index_trace_vir = {};
	if (options_.get_str("CUTOFF_BY") == "THRESHOLD") {
		for (int i = 0; i < noccpi[0]; i++) {
			if (lo->get(0, i) > thresh) {
				index_trace_occ.push_back(i);
				outfile->Printf("\n Occupied orbital %d is partitioned to A with eigenvalue %8.8f", i, lo->get(0, i));
			}
		}
		for (int i = 0; i < nvirpi[0]; i++) {
			if (lv->get(0, i) > thresh) {
				index_trace_vir.push_back(i);
				outfile->Printf("\n Virtual orbital %d is partitioned to A with eigenvalue %8.8f", i, lv->get(0, i));
			}
		}
	}
	
	if (options_.get_str("CUTOFF_BY") == "NUMBER") {
		for (int i = 0; i < num_occ; i++) {
			index_trace_occ.push_back(i);
			outfile->Printf("\n Occupied orbital %d is partitioned to A with eigenvalue %8.8f", i, lo->get(0, i));
		}
		for (int i = 0; i < num_vir; i++) {
			index_trace_vir.push_back(i);
			outfile->Printf("\n Occupied orbital %d is partitioned to A with eigenvalue %8.8f", i, lv->get(0, i));
		}
	}

	//Rotate Bocc to first block, rotate Bvir to last block, in order to freeze them
	//Change to block operation soon
	SharedMatrix Ca_Rt(Ca_t->clone());
	Ca_Rt->zero();

	//Write Bocc
	int sizeBO = noccpi[0] - index_trace_occ.size();
	int sizeBV = nvirpi[0] - index_trace_vir.size();
	int sizeAO = index_trace_occ.size();
	int sizeAV = index_trace_vir.size();
	for (int i = 0; i < sizeBO; ++i) {
		Ca_Rt->set_column(0, i, Ca_->get_column(0, sizeAO + i));
	}

	//Write Aocc
	for (int i = 0; i < sizeAO; ++i) {
		Ca_Rt->set_column(0, i + sizeBO, Ca_->get_column(0, i));
	}

	//Write Avir
	for (int i = 0; i < sizeAV; ++i) {
		Ca_Rt->set_column(0, i + sizeBO + sizeAO, Ca_->get_column(0, i + sizeAO + sizeBO));
	}

	//Write Bvir
	for (int i = 0; i < sizeBV; ++i) {
		Ca_Rt->set_column(0, i + sizeBO + sizeAO + sizeAV, Ca_->get_column(0, sizeBO + sizeAO + sizeAV + i));
	}

	outfile->Printf("\n  The MO coefficients after localization and rotation: \n", sizeBO);
	Ca_Rt->print();
	outfile->Printf("\n  FROZEN_DOCC     = %d", sizeBO);
	outfile->Printf("\n  FROZEN_UOCC	 = %d", sizeBV);

	if (options_.get_bool("Semi-canonicalize A") == true) {
		//Write semi-canonicalization code here
	}

	//Write MO space info and print
	if (options_.get_bool("WRITE_FREEZE_MO") == true) {
		for (int h = 0; h < nirrep_; h++) {
			options_["FROZEN_DOCC"].add(h);
			options_["FROZEN_UOCC"].add(h);
			options_["FROZEN_DOCC"][h].assign(sizeBO);
			options_["FROZEN_UOCC"][h].assign(sizeBV);
		}
	}

	return 0.0;
}

}} // End namespaces

