import os
import glob

def rename_files(directory):
    for filename in glob.glob(directory):
        if os.path.isfile(filename):
            name = filename.split('.fcl')[0]
            ext = '.fcl'
            print(name, ext)
            new_filename = f"{name}_sbnd{ext}"
            new_file_path = os.path.join(directory, new_filename)
            os.system(f'mv {filename} {new_filename}')
            print(f"Renamed: {filename} -> {new_filename}")

rename_files('./*.fcl')
