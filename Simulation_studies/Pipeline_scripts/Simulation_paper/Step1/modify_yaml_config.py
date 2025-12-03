#!/usr/bin/env python3
"""
Helper script to modify YAML configuration file for batch simulations
"""
import sys
import yaml
import os
from pathlib import Path

def modify_yaml_config(template_path, output_path, modifications):
    """
    Modify a YAML config file with given parameters
    
    Args:
        template_path: Path to template YAML file
        output_path: Path to save modified YAML file
        modifications: Dict of parameters to modify
    """
    # Load template config
    with open(template_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Apply modifications
    for key, value in modifications.items():
        config[key] = value
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Save modified config
    with open(output_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    
    print(f"Modified config saved to: {output_path}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python modify_yaml_config.py <template_path> <output_path> <key1>=<value1> [<key2>=<value2> ...]")
        sys.exit(1)
    
    template_path = sys.argv[1]
    output_path = sys.argv[2]
    
    # Parse modifications from command line arguments
    modifications = {}
    for arg in sys.argv[3:]:
        if '=' not in arg:
            continue
        key, value = arg.split('=', 1)
        
        # Try to convert to appropriate type
        try:
            # Try integer first
            value = int(value)
        except ValueError:
            try:
                # Try float
                value = float(value)
            except ValueError:
                # Keep as string
                pass
        
        modifications[key] = value
    
    modify_yaml_config(template_path, output_path, modifications)
