from argparse import ArgumentParser

import ROOT

# /eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_Zbb_ecm91_EvtGen_Bd2KsNuNu
# /eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_Zbb_ecm91_EvtGen_Bd2KstNuNu
# /eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_Zbb_ecm91_EvtGen_Bs2PhiNuNu
# /eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_Zbb_ecm91_EvtGen_Lb2L1520NuNu
# /eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_Zbb_ecm91_EvtGen_Lb2LNuNu

decay_to_candidates = {
    "Bd2KstNuNu": "KPi",
    "Bs2PhiNuNu": "KK",
    "Lb2LNuNu": "pPi",
    "Lb2L1520NuNu": "pK",
}

decay_to_pdgids = {
    "Bd2KstNuNu": ["313", "511"],
    "Bs2PhiNuNu": ["333", "531"],
    "Lb2LNuNu": ["3122", "5122"],
    "Lb2LNuNu": ["3124", "5122"],
}

n_threads = 1
campaign = "winter2023"

if campaign == "spring2021":
    test_files = {
        "Bd2KstNuNu": "events_000143365.root",
        "Bs2PhiNuNu": "events_000027625.root",
        "Lb2LNuNu": "events_000424799.root",
    }
elif campaign == "winter2023":
    test_files = {
        "Bd2KsNuNu": "events_003812794.root",
        "Bd2KstNuNu": "events_006158724.root",
        "Bs2PhiNuNu": "events_005322435.root",
        "Lb2L1520NuNu": "events_011406545.root",
        "Lb2LNuNu": "events_000568214.root"
    }
else:
    raise ValueError(f"Campaign {campaign} has no known test files.")

class Analysis():

    def __init__(self, cmdline_args):
        parser = ArgumentParser(
            description='Additional analysis arguments',
            usage='Provide additional arguments after analysis script path')
        parser.add_argument('--output-dir', default='outputs/stage1', type=str,
                            help='Output directory for the analysis.')
        parser.add_argument('--decay', required=True, choices=list(decay_to_candidates),
                            help='Decay to reconstruct.')
        parser.add_argument('--mva', default='', type=str,
                            help='Path to the trained MVA ROOT file. Without it, tuples are made for BDT training.')
        parser.add_argument('--mva-cut', default=-1., type=float,
                            help='Cut on EVT_MVA1.')
        parser.add_argument('--fraction', default=1., type=float,
                            help='Fraction of each sample to process.')
        parser.add_argument('--chunks', default=1, type=int,
                            help='Number of chunks (batch jobs) per sample.')
        parser.add_argument('--condor', action='store_true', default=False,
                            help='Run on HTCondor.')
        parser.add_argument('--noMVA', default=False,
                            help='Do not run MVA training.')
        self.ana_args, _ = parser.parse_known_args(cmdline_args['unknown'])

        self.decay = self.ana_args.decay
        self.candidates = decay_to_candidates[self.decay]
        self.child_pdgid, self.parent_pdgid = decay_to_pdgids[self.decay]
        self.training = self.ana_args.mva == ''
        self.noMVA = self.ana_args.noMVA 

        # Do only train if training and not noMVA option
        if not self.training and not noMVA:
            ROOT.gInterpreter.ProcessLine(f'''
            TMVA::Experimental::RBDT<> bdt("{self.decay}_BDT", "{self.ana_args.mva}");
            computeModel = TMVA::Experimental::Compute<18, float>(bdt);
            ''')

        sample_opts = {'fraction': self.ana_args.fraction, 'chunks': self.ana_args.chunks}
        self.process_list = {
            f'p8_ee_Zbb_ecm91_EvtGen_{self.decay}': dict(sample_opts),
            'p8_ee_Zbb_ecm91': dict(sample_opts),
            'p8_ee_Zcc_ecm91': dict(sample_opts),
            'p8_ee_Zuds_ecm91': dict(sample_opts),
        }

        self.prod_tag = f'FCCee/{campaign}/IDEA/'

        mode = 'training' if self.training else 'mva'
        self.output_dir = f'{self.ana_args.output_dir}/{mode}/{self.decay}'

        self.run_batch = self.ana_args.condor

        MC_name = 'Lb2LbNuNu' if (self.decay=='Lb2LNuNu' and campaign=='spring2021') else self.decay
        self.test_file = 'root://eospublic.cern.ch//eos/experiment/fcc/ee/' \
                         f'generation/DelphesEvents/{campaign}/IDEA/' \
                         f'p8_ee_Zbb_ecm91_EvtGen_{MC_name}/{test_files[self.decay]}'

    def analyzers(self, dframe):
        decay = self.decay
        candidates = self.candidates
        child_pdgid = self.child_pdgid
        parent_pdgid = self.parent_pdgid

        dframe2 = (dframe
               #############################################
               ##          Aliases for # in python        ##
               #############################################
               .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
               .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
               .Alias("Particle0", "Particle#0.index")
               .Alias("Particle1", "Particle#1.index")


               #############################################
               ##MC record to study the Z->bb events types##
               #############################################               
               .Define("MC_PDG", "FCCAnalyses::MCParticle::get_pdg(Particle)")
               .Define("MC_n",   "int(MC_PDG.size())")
               #.Define("MC_M1",  "FCCAnalyses::myUtils::get_MCMother1(Particle,Particle0)")
               #.Define("MC_M2",  "FCCAnalyses::myUtils::get_MCMother2(Particle,Particle0)")
               #.Define("MC_D1",  "FCCAnalyses::myUtils::get_MCDaughter1(Particle,Particle1)")
               #.Define("MC_D2",  "FCCAnalyses::myUtils::get_MCDaughter2(Particle,Particle1)")
               .Define("MC_M1",  "FCCAnalyses::myUtils::getMC_parent(0,Particle,Particle0)")
               .Define("MC_M2",  "FCCAnalyses::myUtils::getMC_parent(1,Particle,Particle0)")
               .Define("MC_D1",  "FCCAnalyses::myUtils::getMC_daughter(0,Particle,Particle1)")
               .Define("MC_D2",  "FCCAnalyses::myUtils::getMC_daughter(1,Particle,Particle1)")
               .Define("MC_D3",  "FCCAnalyses::myUtils::getMC_daughter(2,Particle,Particle1)")
               .Define("MC_D4",  "FCCAnalyses::myUtils::getMC_daughter(3,Particle,Particle1)")
               .Define("MC_orivtx_x",   "FCCAnalyses::MCParticle::get_vertex_x(Particle)")
               .Define("MC_orivtx_y",   "FCCAnalyses::MCParticle::get_vertex_y(Particle)")
               .Define("MC_orivtx_z",   "FCCAnalyses::MCParticle::get_vertex_z(Particle)")
               .Define("MC_endvtx_x",   "FCCAnalyses::MCParticle::get_endPoint_x(Particle)")
               .Define("MC_endvtx_y",   "FCCAnalyses::MCParticle::get_endPoint_y(Particle)")
               .Define("MC_endvtx_z",   "FCCAnalyses::MCParticle::get_endPoint_z(Particle)")
               .Define("MC_p",   "FCCAnalyses::MCParticle::get_p(Particle)")
               .Define("MC_pt",  "FCCAnalyses::MCParticle::get_pt(Particle)")
               .Define("MC_px",  "FCCAnalyses::MCParticle::get_pt(Particle)")
               .Define("MC_py",  "FCCAnalyses::MCParticle::get_pt(Particle)")
               .Define("MC_pz",  "FCCAnalyses::MCParticle::get_pt(Particle)")
               .Define("MC_e",   "FCCAnalyses::MCParticle::get_e(Particle)")
               .Define("MC_m",   "FCCAnalyses::MCParticle::get_mass(Particle)")
               .Define("MC_q",   "FCCAnalyses::MCParticle::get_charge(Particle)")
               .Define("MC_eta", "FCCAnalyses::MCParticle::get_eta(Particle)")
               .Define("MC_phi", "FCCAnalyses::MCParticle::get_phi(Particle)")
               

               #############################################
               ##               Build MC Vertex           ##
               #############################################
               .Define("MCVertexObject", "FCCAnalyses::myUtils::get_MCVertexObject(Particle, Particle0)")
               .Define("MC_Vertex_x",    "FCCAnalyses::myUtils::get_MCVertex_x(MCVertexObject)")
               .Define("MC_Vertex_y",    "FCCAnalyses::myUtils::get_MCVertex_y(MCVertexObject)")
               .Define("MC_Vertex_z",    "FCCAnalyses::myUtils::get_MCVertex_z(MCVertexObject)")
               .Define("MC_Vertex_ind",  "FCCAnalyses::myUtils::get_MCindMCVertex(MCVertexObject)")
               .Define("MC_Vertex_ntrk", "FCCAnalyses::myUtils::get_NTracksMCVertex(MCVertexObject)")
               .Define("MC_Vertex_n",    "int(MC_Vertex_x.size())")
               .Define("MC_Vertex_PDG",  "FCCAnalyses::myUtils::get_MCpdgMCVertex(MCVertexObject, Particle)")
               .Define("MC_Vertex_PDGmother",  "FCCAnalyses::myUtils::get_MCpdgMotherMCVertex(MCVertexObject, Particle)")
               .Define("MC_Vertex_PDGgmother", "FCCAnalyses::myUtils::get_MCpdgGMotherMCVertex(MCVertexObject, Particle)")
               #############################################
               ##              Build Reco Vertex          ##
               #############################################
               .Define("VertexObject", "FCCAnalyses::myUtils::get_VertexObject(MCVertexObject,ReconstructedParticles,EFlowTrack_1,MCRecoAssociations0,MCRecoAssociations1)")


               #############################################
               ##          Build PV var and filter        ##
               #############################################
               .Define("EVT_hasPV",    "FCCAnalyses::myUtils::hasPV(VertexObject)")
               .Define("EVT_NtracksPV", "float(FCCAnalyses::myUtils::get_PV_ntracks(VertexObject))")
               .Define("EVT_NVertex",   "float(VertexObject.size())")
               .Filter("EVT_hasPV==1")


               #############################################
               ##          Build RECO P with PID          ##
               #############################################
               #.Define("MisIDRate", misid_rate)
               #.Define("NoMisIDPID" ,f"FCCAnalyses::myUtils::PID(ReconstructedParticles, MCRecoAssociations0,MCRecoAssociations1,Particle, 0.)")
               #.Define("RecoPartPID" ,f"FCCAnalyses::myUtils::PID(ReconstructedParticles, MCRecoAssociations0,MCRecoAssociations1,Particle, {misid_rate})")
               .Define("RecoPartPID" ,"FCCAnalyses::myUtils::PID(ReconstructedParticles, MCRecoAssociations0,MCRecoAssociations1,Particle)")
               

               #############################################
               ##    Build RECO P with PID at vertex      ##
               #############################################
               .Define("RecoPartPIDAtVertex" ,"FCCAnalyses::myUtils::get_RP_atVertex(RecoPartPID, VertexObject)")


               #############################################
               ##         Build vertex variables          ##
               #############################################
               .Define("Vertex_x",        "FCCAnalyses::myUtils::get_Vertex_x(VertexObject)")
               .Define("Vertex_y",        "FCCAnalyses::myUtils::get_Vertex_y(VertexObject)")
               .Define("Vertex_z",        "FCCAnalyses::myUtils::get_Vertex_z(VertexObject)")
               .Define("Vertex_xErr",     "FCCAnalyses::myUtils::get_Vertex_xErr(VertexObject)")
               .Define("Vertex_yErr",     "FCCAnalyses::myUtils::get_Vertex_yErr(VertexObject)")
               .Define("Vertex_zErr",     "FCCAnalyses::myUtils::get_Vertex_zErr(VertexObject)")

               .Define("Vertex_chi2",     "FCCAnalyses::myUtils::get_Vertex_chi2(VertexObject)")
               .Define("Vertex_mcind",    "FCCAnalyses::myUtils::get_Vertex_indMC(VertexObject)")
               .Define("Vertex_ind",      "FCCAnalyses::myUtils::get_Vertex_ind(VertexObject)")
               .Define("Vertex_isPV",     "FCCAnalyses::myUtils::get_Vertex_isPV(VertexObject)")
               .Define("Vertex_ntrk",     "FCCAnalyses::myUtils::get_Vertex_ntracks(VertexObject)")
               .Define("Vertex_n",        "int(Vertex_x.size())")
               .Define("Vertex_mass",     "FCCAnalyses::myUtils::get_Vertex_mass(VertexObject,RecoPartPIDAtVertex)")

               .Define("Vertex_d2PV",     "FCCAnalyses::myUtils::get_Vertex_d2PV(VertexObject,-1)")
               .Define("Vertex_d2PVx",    "FCCAnalyses::myUtils::get_Vertex_d2PV(VertexObject,0)")
               .Define("Vertex_d2PVy",    "FCCAnalyses::myUtils::get_Vertex_d2PV(VertexObject,1)")
               .Define("Vertex_d2PVz",    "FCCAnalyses::myUtils::get_Vertex_d2PV(VertexObject,2)")
               
               .Define("Vertex_d2PVErr",  "FCCAnalyses::myUtils::get_Vertex_d2PVError(VertexObject,-1)")
               .Define("Vertex_d2PVxErr", "FCCAnalyses::myUtils::get_Vertex_d2PVError(VertexObject,0)")
               .Define("Vertex_d2PVyErr", "FCCAnalyses::myUtils::get_Vertex_d2PVError(VertexObject,1)")
               .Define("Vertex_d2PVzErr", "FCCAnalyses::myUtils::get_Vertex_d2PVError(VertexObject,2)")
               
               .Define("Vertex_d2PVSig",  "Vertex_d2PV/Vertex_d2PVErr")
               .Define("Vertex_d2PVxSig", "Vertex_d2PVx/Vertex_d2PVxErr")
               .Define("Vertex_d2PVySig", "Vertex_d2PVy/Vertex_d2PVyErr")
               .Define("Vertex_d2PVzSig", "Vertex_d2PVz/Vertex_d2PVzErr")

               .Define("Vertex_d2MC",     "FCCAnalyses::myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,-1)")
               .Define("Vertex_d2MCx",    "FCCAnalyses::myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,0)")
               .Define("Vertex_d2MCy",    "FCCAnalyses::myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,1)")
               .Define("Vertex_d2MCz",    "FCCAnalyses::myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,2)")

               .Define("EVT_dPV2DVmin",   "FCCAnalyses::myUtils::get_dPV2DV_min(Vertex_d2PV)")
               .Define("EVT_dPV2DVmax",   "FCCAnalyses::myUtils::get_dPV2DV_max(Vertex_d2PV)")
               .Define("EVT_dPV2DVave",   "FCCAnalyses::myUtils::get_dPV2DV_ave(Vertex_d2PV)")
               

               #############################################
               ##        Build Kstz -> KPi  candidates      ##
               #############################################
               .Define(f"{candidates}Candidates",         f"FCCAnalyses::myUtils::build_{decay}(VertexObject,RecoPartPIDAtVertex)")


               #############################################
               ##       Filter Kstz -> KPi candidates      ##
               ############################################# 
               .Define(f"EVT_N{candidates}",              f"float(FCCAnalyses::myUtils::getFCCAnalysesComposite_N({candidates}Candidates))")
               .Filter(f"EVT_N{candidates}>0")


               #############################################
               ##    Attempt to add a truth match         ##
               #############################################
               #.Define("TruthMatching" ,f"FCCAnalyses::myUtils::add_truthmatched2({candidates}Candidates, MCParticles, MCRecoAssociations0, ReconstructedParticles, MCRecoAssociations1)")
               .Define("TruthMatching" ,f"FCCAnalyses::myUtils::add_truthmatched2({candidates}Candidates, Particle, VertexObject, MCRecoAssociations1, ReconstructedParticles, Particle0)")


               #############################################
               ##              Build the thrust           ##
               ############################################# 
               .Define("RP_e",          "FCCAnalyses::ReconstructedParticle::get_e(RecoPartPIDAtVertex)")
               .Define("RP_px",         "FCCAnalyses::ReconstructedParticle::get_px(RecoPartPIDAtVertex)")
               .Define("RP_py",         "FCCAnalyses::ReconstructedParticle::get_py(RecoPartPIDAtVertex)")
               .Define("RP_pz",         "FCCAnalyses::ReconstructedParticle::get_pz(RecoPartPIDAtVertex)")
               .Define("RP_charge",     "FCCAnalyses::ReconstructedParticle::get_charge(RecoPartPIDAtVertex)")
              
               .Define("EVT_thrustNP",      'FCCAnalyses::Algorithms::minimize_thrust("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
               .Define("RP_thrustangleNP",  'FCCAnalyses::Algorithms::getAxisCosTheta(EVT_thrustNP, RP_px, RP_py, RP_pz)')
               .Define("EVT_thrust",        'FCCAnalyses::Algorithms::getThrustPointing(1.)(RP_thrustangleNP, RP_e, EVT_thrustNP)') # changed from 'Algorithms::getThrustPointing(RP_thrustangleNP, RP_e, EVT_thrustNP, 1.)' because of https://github.com/HEP-FCC/FCCAnalyses/commit/e9c4787f82505115be0c084da4453031c3cf8fdf
               .Define("RP_thrustangle",    'FCCAnalyses::Algorithms::getAxisCosTheta(EVT_thrust, RP_px, RP_py, RP_pz)')

               
               #############################################
               ##        Get thrust related values        ##
               ############################################# 
               ##hemis0 == negative angle == max energy hemisphere if pointing
               ##hemis1 == positive angle == min energy hemisphere if pointing
               .Define("EVT_thrusthemis0_n",    "FCCAnalyses::Algorithms::getAxisN(0)(RP_thrustangle, RP_charge)")
               .Define("EVT_thrusthemis1_n",    "FCCAnalyses::Algorithms::getAxisN(1)(RP_thrustangle, RP_charge)")
               .Define("EVT_thrusthemis0_e",    "FCCAnalyses::Algorithms::getAxisEnergy(0)(RP_thrustangle, RP_charge, RP_e)")
               .Define("EVT_thrusthemis1_e",    "FCCAnalyses::Algorithms::getAxisEnergy(1)(RP_thrustangle, RP_charge, RP_e)")

               .Define("EVT_ThrustEmax_E",         "EVT_thrusthemis0_e.at(0)")
               .Define("EVT_ThrustEmax_Echarged",  "EVT_thrusthemis0_e.at(1)")
               .Define("EVT_ThrustEmax_Eneutral",  "EVT_thrusthemis0_e.at(2)")
               .Define("EVT_ThrustEmax_N",         "float(EVT_thrusthemis0_n.at(0))")
               .Define("EVT_ThrustEmax_Ncharged",  "float(EVT_thrusthemis0_n.at(1))")
               .Define("EVT_ThrustEmax_Nneutral",  "float(EVT_thrusthemis0_n.at(2))")

               .Define("EVT_ThrustEmin_E",         "EVT_thrusthemis1_e.at(0)")
               .Define("EVT_ThrustEmin_Echarged",  "EVT_thrusthemis1_e.at(1)")
               .Define("EVT_ThrustEmin_Eneutral",  "EVT_thrusthemis1_e.at(2)")
               .Define("EVT_ThrustEmin_N",         "float(EVT_thrusthemis1_n.at(0))")
               .Define("EVT_ThrustEmin_Ncharged",  "float(EVT_thrusthemis1_n.at(1))")
               .Define("EVT_ThrustEmin_Nneutral",  "float(EVT_thrusthemis1_n.at(2))")


               .Define("Vertex_thrust_angle",   "FCCAnalyses::myUtils::get_Vertex_thrusthemis_angle(VertexObject, RecoPartPIDAtVertex, EVT_thrust)")
               .Define("DVertex_thrust_angle",  "FCCAnalyses::myUtils::get_DVertex_thrusthemis_angle(VertexObject, RecoPartPIDAtVertex, EVT_thrust)")
               ###0 == negative angle==max energy , 1 == positive angle == min energy
               .Define("Vertex_thrusthemis_emin",    "FCCAnalyses::myUtils::get_Vertex_thrusthemis(Vertex_thrust_angle, 1)")
               .Define("Vertex_thrusthemis_emax",    "FCCAnalyses::myUtils::get_Vertex_thrusthemis(Vertex_thrust_angle, 0)")

               .Define("EVT_ThrustEmin_NDV", "float(FCCAnalyses::myUtils::get_Npos(DVertex_thrust_angle))")
               .Define("EVT_ThrustEmax_NDV", "float(FCCAnalyses::myUtils::get_Nneg(DVertex_thrust_angle))")

               .Define("EVT_Thrust_Mag",  "EVT_thrust.at(0)")
               .Define("EVT_Thrust_X",    "EVT_thrust.at(1)")
               .Define("EVT_Thrust_XErr", "EVT_thrust.at(2)")
               .Define("EVT_Thrust_Y",    "EVT_thrust.at(3)")
               .Define("EVT_Thrust_YErr", "EVT_thrust.at(4)")
               .Define("EVT_Thrust_Z",    "EVT_thrust.at(5)")
               .Define("EVT_Thrust_ZErr", "EVT_thrust.at(6)")


               .Define("DV_tracks", "FCCAnalyses::myUtils::get_pseudotrack(VertexObject,RecoPartPIDAtVertex)")

               .Define("DV_d0",            "FCCAnalyses::myUtils::get_trackd0(DV_tracks)")
               .Define("DV_z0",            "FCCAnalyses::myUtils::get_trackz0(DV_tracks)")

               .Define(f"{candidates}Candidates_mass",    f"FCCAnalyses::myUtils::getFCCAnalysesComposite_mass({candidates}Candidates)")
               .Define(f"{candidates}Candidates_q",       f"FCCAnalyses::myUtils::getFCCAnalysesComposite_charge({candidates}Candidates)")
               .Define(f"{candidates}Candidates_vertex",  f"FCCAnalyses::myUtils::getFCCAnalysesComposite_vertex({candidates}Candidates)")
               .Define(f"{candidates}Candidates_mcvertex",f"FCCAnalyses::myUtils::getFCCAnalysesComposite_mcvertex({candidates}Candidates,VertexObject)")
               .Define(f"{candidates}Candidates_truth",   f"FCCAnalyses::myUtils::getFCCAnalysesComposite_truthMatch(TruthMatching)") # KPiCandidates -> TruthMatching
               .Define(f"{candidates}Candidates_px",      f"FCCAnalyses::myUtils::getFCCAnalysesComposite_p({candidates}Candidates,0)")
               .Define(f"{candidates}Candidates_py",      f"FCCAnalyses::myUtils::getFCCAnalysesComposite_p({candidates}Candidates,1)")
               .Define(f"{candidates}Candidates_pz",      f"FCCAnalyses::myUtils::getFCCAnalysesComposite_p({candidates}Candidates,2)")
               .Define(f"{candidates}Candidates_p",       f"FCCAnalyses::myUtils::getFCCAnalysesComposite_p({candidates}Candidates,-1)")
               .Define(f"{candidates}Candidates_B",       f"FCCAnalyses::myUtils::getFCCAnalysesComposite_B({candidates}Candidates, VertexObject, RecoPartPIDAtVertex)")
               
               .Define(f"{candidates}Candidates_track",   f"FCCAnalyses::myUtils::getFCCAnalysesComposite_track({candidates}Candidates, VertexObject)")
               .Define(f"{candidates}Candidates_d0",      f"FCCAnalyses::myUtils::get_trackd0({candidates}Candidates_track)")
               .Define(f"{candidates}Candidates_z0",      f"FCCAnalyses::myUtils::get_trackz0({candidates}Candidates_track)")

               .Define(f"{candidates}Candidates_anglethrust", f"FCCAnalyses::myUtils::getFCCAnalysesComposite_anglethrust({candidates}Candidates, EVT_thrust)")
               .Define("CUT_hasCandEmin",           f"FCCAnalyses::myUtils::has_anglethrust_emin({candidates}Candidates_anglethrust)")
               .Filter("CUT_hasCandEmin>0")
               
               .Define(f"{candidates}Candidates_h1px",   f"FCCAnalyses::myUtils::getFCCAnalysesComposite_p({candidates}Candidates, VertexObject, RecoPartPIDAtVertex, 0, 0)")
               .Define(f"{candidates}Candidates_h1py",   f"FCCAnalyses::myUtils::getFCCAnalysesComposite_p({candidates}Candidates, VertexObject, RecoPartPIDAtVertex, 0, 1)")
               .Define(f"{candidates}Candidates_h1pz",   f"FCCAnalyses::myUtils::getFCCAnalysesComposite_p({candidates}Candidates, VertexObject, RecoPartPIDAtVertex, 0, 2)")
               .Define(f"{candidates}Candidates_h1p",    f"FCCAnalyses::myUtils::getFCCAnalysesComposite_p({candidates}Candidates, VertexObject, RecoPartPIDAtVertex, 0, -1)")
               .Define(f"{candidates}Candidates_h1q",    f"FCCAnalyses::myUtils::getFCCAnalysesComposite_q({candidates}Candidates, VertexObject, RecoPartPIDAtVertex, 0)")
               .Define(f"{candidates}Candidates_h1m",    f"FCCAnalyses::myUtils::getFCCAnalysesComposite_mass({candidates}Candidates, VertexObject, RecoPartPIDAtVertex, 0)")
               .Define(f"{candidates}Candidates_h1type", f"FCCAnalyses::myUtils::getFCCAnalysesComposite_type({candidates}Candidates, VertexObject, RecoPartPIDAtVertex, 0)")
               .Define(f"{candidates}Candidates_h1d0",   f"FCCAnalyses::myUtils::getFCCAnalysesComposite_d0({candidates}Candidates, VertexObject, 0)")
               .Define(f"{candidates}Candidates_h1z0",   f"FCCAnalyses::myUtils::getFCCAnalysesComposite_z0({candidates}Candidates, VertexObject, 0)")
               
               .Define(f"{candidates}Candidates_h2px",   f"FCCAnalyses::myUtils::getFCCAnalysesComposite_p({candidates}Candidates, VertexObject, RecoPartPIDAtVertex, 1, 0)")
               .Define(f"{candidates}Candidates_h2py",   f"FCCAnalyses::myUtils::getFCCAnalysesComposite_p({candidates}Candidates, VertexObject, RecoPartPIDAtVertex, 1, 1)")
               .Define(f"{candidates}Candidates_h2pz",   f"FCCAnalyses::myUtils::getFCCAnalysesComposite_p({candidates}Candidates, VertexObject, RecoPartPIDAtVertex, 1, 2)")
               .Define(f"{candidates}Candidates_h2p",    f"FCCAnalyses::myUtils::getFCCAnalysesComposite_p({candidates}Candidates, VertexObject, RecoPartPIDAtVertex, 1, -1)")
               .Define(f"{candidates}Candidates_h2q",    f"FCCAnalyses::myUtils::getFCCAnalysesComposite_q({candidates}Candidates, VertexObject, RecoPartPIDAtVertex, 1)")
               .Define(f"{candidates}Candidates_h2m",    f"FCCAnalyses::myUtils::getFCCAnalysesComposite_mass({candidates}Candidates, VertexObject, RecoPartPIDAtVertex, 1)")
               .Define(f"{candidates}Candidates_h2type", f"FCCAnalyses::myUtils::getFCCAnalysesComposite_type({candidates}Candidates, VertexObject, RecoPartPIDAtVertex, 1)")
               .Define(f"{candidates}Candidates_h2d0",   f"FCCAnalyses::myUtils::getFCCAnalysesComposite_d0({candidates}Candidates, VertexObject, 1)")
               .Define(f"{candidates}Candidates_h2z0",   f"FCCAnalyses::myUtils::getFCCAnalysesComposite_z0({candidates}Candidates, VertexObject, 1)")
               
               .Define(f"True{candidates}_vertex",        f"FCCAnalyses::myUtils::get_trueVertex(MCVertexObject,Particle,Particle0, {child_pdgid}, {parent_pdgid})")
               .Define(f"True{candidates}_track",         f"FCCAnalyses::myUtils::get_truetrack(True{candidates}_vertex, MCVertexObject, Particle)")
               .Define(f"True{candidates}_d0",            f"FCCAnalyses::myUtils::get_trackd0(True{candidates}_track)")
               .Define(f"True{candidates}_z0",            f"FCCAnalyses::myUtils::get_trackz0(True{candidates}_track)")
           )

        if not self.training and not noMVA:
            dframe2 = (dframe2
               # Build MVA 
               .Define("MVAVec", ROOT.computeModel, ("EVT_ThrustEmin_E",        "EVT_ThrustEmax_E",
                                                     "EVT_ThrustEmin_Echarged", "EVT_ThrustEmax_Echarged",
                                                     "EVT_ThrustEmin_Eneutral", "EVT_ThrustEmax_Eneutral",
                                                     "EVT_ThrustEmin_Ncharged", "EVT_ThrustEmax_Ncharged",
                                                     "EVT_ThrustEmin_Nneutral", "EVT_ThrustEmax_Nneutral",
                                                     "EVT_NtracksPV",           "EVT_NVertex",
                                                     f"EVT_N{candidates}",                "EVT_ThrustEmin_NDV",
                                                     "EVT_ThrustEmax_NDV",      "EVT_dPV2DVmin",
                                                     "EVT_dPV2DVmax",           "EVT_dPV2DVave"))
               .Define("EVT_MVA1", "MVAVec.at(0)")
               .Filter(f"EVT_MVA1>{self.ana_args.mva_cut}")
            )

        return dframe2

    ## Felicia: Check if all needed variables in here. 
    def output(self):
        candidates = self.candidates
        branch_list = [
                "MC_PDG","MC_M1","MC_M2","MC_n","MC_D1","MC_D2","MC_D3","MC_D4",
                "MC_p","MC_pt","MC_px","MC_py","MC_pz","MC_eta","MC_phi",
                "MC_orivtx_x","MC_orivtx_y","MC_orivtx_z", 
                "MC_endvtx_x", "MC_endvtx_y", "MC_endvtx_z", "MC_e","MC_m",
                "EVT_ThrustEmin_E",          "EVT_ThrustEmax_E",
                "EVT_ThrustEmin_Echarged",   "EVT_ThrustEmax_Echarged",
                "EVT_ThrustEmin_Eneutral",   "EVT_ThrustEmax_Eneutral",
                "EVT_ThrustEmin_N",          "EVT_ThrustEmax_N",                
                "EVT_ThrustEmin_Ncharged",   "EVT_ThrustEmax_Ncharged",
                "EVT_ThrustEmin_Nneutral",   "EVT_ThrustEmax_Nneutral",
                "EVT_ThrustEmin_NDV",        "EVT_ThrustEmax_NDV",
                "EVT_Thrust_Mag",
                "EVT_Thrust_X",  "EVT_Thrust_XErr",
                "EVT_Thrust_Y",  "EVT_Thrust_YErr",
                "EVT_Thrust_Z",  "EVT_Thrust_ZErr",

                "EVT_NtracksPV", "EVT_NVertex", f"EVT_N{candidates}",
                
                "EVT_dPV2DVmin","EVT_dPV2DVmax","EVT_dPV2DVave",

                "MC_Vertex_x", "MC_Vertex_y", "MC_Vertex_z", 
                "MC_Vertex_ntrk", "MC_Vertex_n",
                
                "MC_Vertex_PDG","MC_Vertex_PDGmother","MC_Vertex_PDGgmother",
                
                "Vertex_x", "Vertex_y", "Vertex_z",
                "Vertex_xErr", "Vertex_yErr", "Vertex_zErr",
                "Vertex_isPV", "Vertex_ntrk", "Vertex_chi2", "Vertex_n",
                "Vertex_thrust_angle", "Vertex_thrusthemis_emin", "Vertex_thrusthemis_emax",

                "Vertex_d2PV", "Vertex_d2PVx", "Vertex_d2PVy", "Vertex_d2PVz",
                "Vertex_d2PVErr", "Vertex_d2PVxErr", "Vertex_d2PVyErr", "Vertex_d2PVzErr",
                "Vertex_mass",
                "DV_d0","DV_z0",
                
                f"True{candidates}_vertex", f"True{candidates}_d0", f"True{candidates}_z0", 
                
                f"{candidates}Candidates_mass", f"{candidates}Candidates_vertex", f"{candidates}Candidates_mcvertex", f"{candidates}Candidates_B",
                f"{candidates}Candidates_truth",
                f"{candidates}Candidates_px", f"{candidates}Candidates_py", f"{candidates}Candidates_pz", f"{candidates}Candidates_p", f"{candidates}Candidates_q",
                f"{candidates}Candidates_d0",  f"{candidates}Candidates_z0",f"{candidates}Candidates_anglethrust",
                
                f"{candidates}Candidates_h1px", f"{candidates}Candidates_h1py", f"{candidates}Candidates_h1pz",
                f"{candidates}Candidates_h1p", f"{candidates}Candidates_h1q", f"{candidates}Candidates_h1m", f"{candidates}Candidates_h1type",
                f"{candidates}Candidates_h1d0", f"{candidates}Candidates_h1z0",
                f"{candidates}Candidates_h2px", f"{candidates}Candidates_h2py", f"{candidates}Candidates_h2pz",
                f"{candidates}Candidates_h2p", f"{candidates}Candidates_h2q", f"{candidates}Candidates_h2m", f"{candidates}Candidates_h2type",
                f"{candidates}Candidates_h2d0", f"{candidates}Candidates_h2z0",
                ]
        if not self.training:
            branch_list.append("EVT_MVA1")
        return branch_list
