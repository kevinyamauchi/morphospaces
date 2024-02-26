import os
import json

import numpy as np
import mrcfile

from tqdm import tqdm
import shutil

# directory to save all of the data to
base_directory = "./data"

# Define the annotation directory
annotations_directory = os.path.join(base_directory, "annotations")
tomo_directory = os.path.join(base_directory, "images")

labels = {
    "ribosome": {"idx": 20, "radius": 6},
    "fatty_acid_synthase": {"idx": 21, "radius": 8},
}


def insert_spheres(data, coordinates, sphere_radius, sphere_value):
    """
    Insert spheres into a volume.
    """
    # Create a meshgrid of coordinates
    x = np.arange(data.shape[0])
    y = np.arange(data.shape[1])
    z = np.arange(data.shape[2])

    # Iterate over the coordinates
    for x, y, z in tqdm(coordinates):
        x,y,z = int(x), int(y), int(z)

        for i in range(-sphere_radius, sphere_radius):
            for j in range(-sphere_radius, sphere_radius):
                for k in range(-sphere_radius, sphere_radius):
                    if (x+i) < 0 or (y+j) < 0 or (z+k) < 0 or (x+i) >= data.shape[0] or (y+j) >= data.shape[1] or (z+k) >= data.shape[2]:
                        continue
                    if np.sqrt(i**2 + j**2 + k**2) < sphere_radius:
                        data[x+i, y+j, z+k] = sphere_value


    return data


## Assumes folder structure:
# data
# ├── annotations
# │   ├── dataset1
# │   │   ├── annotation1.ndjson
# │   │   ├── annotation2.ndjson
# │   │   └── segmentation_mask1.mrc
# │   │   └── ...
# │   └── dataset2
# │       ├── annotation1.ndjson
# │       ├── annotation2.ndjson
# │       └── segmentation_mask1.mrc
# │       └── ...
# └── images
#     ├── dataset1.mrc
#     └── dataset2.mrc


# Iterate over the datasets (each dataset is a tomogram with associated annotations)
for dataset_name in os.listdir(annotations_directory):
    if not os.path.isdir(os.path.join(annotations_directory, dataset_name)):
        continue
    tomo_annotations_directory = os.path.join(
        annotations_directory, dataset_name
    )

    # Load the first annotation (all annotation files contain the entire segmentation data)
    for filename in os.listdir(tomo_annotations_directory):
        if filename.endswith(".mrc") and not filename == "merged_point_annotations.mrc":
            # Load the annotation
            with mrcfile.open(os.path.join(tomo_annotations_directory, filename), mode="r", permissive=True) as mrc:
                # load the annotation data
                annotation = mrc.data.copy()

                # Transpose the data to match the coordinates
                annotation = np.transpose(annotation, (2, 1, 0))
            break

    # Iterate over the annotation files
    for filename in os.listdir(tomo_annotations_directory):

        if filename.endswith(".ndjson"):
            # Extract the protein name from the filename
            for protein in labels.keys():
                if protein in filename:
                    break
            
            # Load the point annotations
            coordinates = []
            with open(os.path.join(tomo_annotations_directory, filename), "r") as f:
                for line in f:
                    # Parse the JSON object from each line
                    data = json.loads(line)
                    
                    # Extract the location object
                    location = data.get('location', {})
                    
                    # Extract x, y, z coordinates
                    x = location.get('x')
                    y = location.get('y')
                    z = location.get('z')
                    
                    # # Append the coordinates to the list
                    coordinates.append((x, y, z))
            
            # Convert the coordinates to a numpy array
            coordinates = np.array(coordinates)

            # Insert the spheres into the annotation volume
            annotation = insert_spheres(
                data=annotation, 
                coordinates=coordinates, 
                sphere_radius=labels[protein]["radius"], 
                sphere_value=labels[protein]["idx"]
                )
            
    # Save the annotation
    output_filename = os.path.join(tomo_annotations_directory, "merged_point_annotations.mrc")
    annotation = np.transpose(annotation, (2, 1, 0))

    # Copying and changing the annotation keeps the original header information
    shutil.copyfile(os.path.join(tomo_annotations_directory, filename), output_filename)
    # change data to the new annotation
    with mrcfile.open(output_filename, mode="r+") as mrc:
        mrc.set_data(annotation)
        mrc.close()
