import os
import json
import numpy as np
import argparse

def load_decoder(N, differentiate, method, dir_WAs):
    """Load precomputed decoder dictionary from JSON file"""
    # Construct file path
    decoder_path = os.path.join(
        dir_WAs,
        f"N_{N}",
        f"diff_{differentiate}",
        f"decoder_{method}.json"
    )
    
    # Validate file existence
    if not os.path.exists(decoder_path):
        raise FileNotFoundError(f"Decoder file missing: {decoder_path}")
    
    # Load and return decoder
    with open(decoder_path, 'r') as f:
        return json.load(f)

def decode_experiment(readout, decoder_dict):
    """Decode experimental readout using precomputed decoder"""
    # Convert readout to tuple key
    readout_tuple = tuple(map(int, readout.astype(bool).tolist()))
    
    # Validate key format
    if not all(isinstance(x, int) for x in readout_tuple):
        raise ValueError("Readout must be convertible to tuple of integers")
    
    # Perform lookup
    if readout_tuple not in decoder_dict:
        # Try alternative representation if needed
        alt_tuple = tuple(1 if x else 0 for x in readout_tuple)
        if alt_tuple in decoder_dict:
            return decoder_dict[alt_tuple]
        raise KeyError(f"Readout pattern {readout_tuple} not found in decoder")
    
    return decoder_dict[readout_tuple]

def main():
    parser = argparse.ArgumentParser(description='Decode pooling experiment results')
    parser.add_argument('--readout', type=str, required=True,
                        help='Comma-separated readout values (e.g., "1,0,1,1")')
    parser.add_argument('--N', type=int, required=True,
                        help='Number of compounds')
    parser.add_argument('--differentiate', type=int, required=True,
                        help='Differentiation level used')
    parser.add_argument('--method', type=str, required=True,
                        help='Assignment method (e.g., random, grid)')
    parser.add_argument('--dir_WAs', type=str, required=True,
                        help='Directory containing decoder files')
    
    args = parser.parse_args()
    
    # Process readout input
    try:
        readout_array = np.array([int(x) for x in args.readout.split(',')])
    except ValueError:
        raise ValueError("Readout must be comma-separated integers (e.g., '1,0,1,1')")
    
    # Load decoder dictionary
    decoder = load_decoder(
        N=args.N,
        differentiate=args.differentiate,
        method=args.method,
        dir_WAs=args.dir_WAs
    )
    
    # Perform decoding
    positives = decode_experiment(readout_array, decoder)
    
    # Format output
    print("\n🔍 Decoding Results:")
    print(f"• Compounds: {args.N}, Differentiation: {args.differentiate}")
    print(f"• Method: {args.method}, Readout: {args.readout}")
    print(f"• Positive compounds: {positives}")

if __name__ == "__main__":
    main()
