import os
import argparse
import pyvista as pv
import numpy as np


# USE: 
# python3 automate_extraction.py "read_folder_path" "write_folder_path"


# Function to extract points and quantities along the line segment and write to file
def extract_points_along_line(file_path, point1, point2, output_folder_vtp, output_folder_csv, sample_points):
    data = pv.read(file_path)

    # Create a line and sample points along it
    sampled_points = data.sample_over_line(point1, point2, sample_points)

    # Create a line and sample points along it
    #sampled_points = data.sample_over_line(point1, point2)
    
    # Convert sampled points to numpy array
    sampled_points_array = sampled_points.points
    
    # Extract scalar/vector field data
    field_data = [sampled_points.point_data[key] for key in sampled_points.point_data.keys()]
    field_data_array = np.column_stack([sampled_points_array] + field_data)
    
    header = 'x,y,z,' + ','.join(sampled_points.point_data.keys())
    
    # Define the output file paths
    output_file_path_vtp = os.path.join(output_folder_vtp, os.path.basename(file_path).replace('.pvti', '_line.vtp'))
    output_file_path_csv = os.path.join(output_folder_csv, os.path.basename(file_path).replace('.pvti', '_line.csv'))
    
    # Save the sampled points to a new .vtp file
    sampled_points.save(output_file_path_vtp)
    print(f"Extracted line segment to: {output_file_path_vtp}")
    
    # Save the sampled points array to a CSV file
    np.savetxt(output_file_path_csv, field_data_array, delimiter=',', header=header, comments='')
    print(f"Extracted field data to:   {output_file_path_csv}")


    print("")




def main():
    # Define the folder containing the .pvti files
    # python3 automate_extraction.py "read_folder_path" "write_folder_path"
    parser = argparse.ArgumentParser(description='Process .pvti files and extract points along a line.')
    parser.add_argument('read_folder_path', type=str, help='Folder containing the .pvti files')
    parser.add_argument('write_folder_path', type=str, help='Folder to write the output files')
    parser.add_argument('sample_points', type=int, help='Number of sample points in data extraction')

    parser.add_argument('point1', type=float, nargs=3, help='Coordinates of point 1 as three float values x,y,z')
    parser.add_argument('point2', type=float, nargs=3, help='Coordinates of point 2 as three float values x,y,z')

    args = parser.parse_args()

    read_folder_path = args.read_folder_path
    write_folder_path = args.write_folder_path
    sample_points = args.sample_points

    point1 = np.array(args.point1)
    point2 = np.array(args.point2)

    # Create output folders if they don't exist
    output_folder_vtp = os.path.join(write_folder_path, 'extracted_lines_vtp')
    os.makedirs(output_folder_vtp, exist_ok=True)

    output_folder_csv = os.path.join(write_folder_path, 'extracted_lines_csv')
    os.makedirs(output_folder_csv, exist_ok=True)

    # Iterate over all .pvti files in the folder and extract points along the line segment
    for file_name in os.listdir(read_folder_path):
        if file_name.endswith('.pvti'):
            read_file_path  = os.path.join( read_folder_path, file_name)
            extract_points_along_line(read_file_path, point1, point2, output_folder_vtp, output_folder_csv, sample_points)


if __name__=="__main__":
    main()
