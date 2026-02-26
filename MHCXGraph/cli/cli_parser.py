import argparse

def none_or_float(value):
    if value == 'None':
        return None
    return float(value)

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', '0'):
        return False
    raise argparse.ArgumentTypeError('Boolean value expected.')

def parse_args():
    """Configura e retorna os argumentos CLI."""
    parser = argparse.ArgumentParser(description="Building common subgraphs")
    parser.add_argument("--manifest", type=str, required=True,
                help="Path to the unified JSON manifest.")
    # parser.add_argument('--folder_path', type=str, default=None,
    #                     help='Folder path with PDB input files.')
    # parser.add_argument('--files_name', type=str, default=None,
    #                     help="PDB files separated by comma ','. If not provided, the system will ask for user to choose in menu")
    # parser.add_argument('--reference_graph', type=str, default=None,
    #                     help="Reference Graph to be used.")
    # parser.add_argument('--interface_list', type=str, default='',
    #                     help='File with a canonical list of MHC residues at the interface with TCR. No header needed for this file.')
    # parser.add_argument('--run_name', type=str, default='test',
    #                     help='Name for storing results in the output folder')
    # parser.add_argument('--association_mode', type=str, default='identity',
    #                     help='Mode for creating association nodes. Identify or similarity.')  
    # parser.add_argument('--exclude_waters', type=str2bool, default=True,
    #                     help='Exclude water molecules from the graphs.')    
    # parser.add_argument('--track_residues', type=str, default=None,
    #                     help='List of residues to track in the process of association. If not provided, none will be tracked.')                                  
    # parser.add_argument('--track_steps', type=str2bool, default=False,
    #                     help='Track steps in the association process. If True, it will track the steps, otherwise it will not.')
    # parser.add_argument('--centroid_granularity', type=str, default='all_atoms',
    #                     help='Granularity for building the graphs. Options are: all_atoms, backbone, side_chain, ca_only.')
    
    # parser.add_argument('--output_path', type=str, default='~/',
    #                     help='Path to store output results.')
    # parser.add_argument('--factors_path', type=str, default=None,
    #                     help="Factors for calculating the residue similarity ")
    # parser.add_argument('--classes_path', type=str, default=None,
    #                     help="JSON File with the residues, RSA or depth agrouped by classes.")
    
    # parser.add_argument('--manifest', type=str, default=None,
    #                     help="Path to Json file which contains the list of constrains to be used for each pdb file")

    # parser.add_argument('--rsa_filter', type=none_or_float, default=0.1,
    #                     help="Threshold for filter residues by RSA")
    # parser.add_argument('--depth_filter', type=none_or_float, default=10,
    #                     help="Threshold for filter residues by Depth")

    # parser.add_argument('--distance_diff_threshold', type=float, default=2.0,
    #                     help="Threshold for distance difference")
    # parser.add_argument('--centroid_threshold', type=float, default=8.5,
    #                     help="Distance threshold for building the molA and molB interface graphs")

    # parser.add_argument('--rsa_bins', type=int, default=5,
    #                     help="The number of bins that the RSA will be equally dicretized")
    # parser.add_argument('--depth_bins', type=int, default=5,
    #                     help="The number of bins that the depth will be equally discretized")
    # parser.add_argument('--distance_bins', type=int, default=5,
    #                     help="The number of bins that the distance will be equally discretized")
 
    # parser.add_argument('--serd_config', type=str, default=None,
    #                     help="Path to Json file which contains the SERD configuration" )
    # parser.add_argument('--angle_diff', type=float, default=20.0,
    #                     help="Max angle difference to filter association graph's nodes")

    # parser.add_argument('--check_depth', type=str2bool, default=True,
    #                     help="Check depth similarity")
    # parser.add_argument('--check_rsa', type=str2bool, default=True,
    #                     help="Check rsa similarity")

    # parser.add_argument('--debug', type=str2bool, default=False,
    #                     help="Activate debug mode")
    return parser.parse_args()
