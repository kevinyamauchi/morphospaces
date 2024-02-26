import os

from cryoet_data_portal import Client, Run
from tqdm import tqdm

# directory to save all of the data to
base_directory = "./data"

# make the output directories
image_directory = os.path.join(base_directory, "images")
os.makedirs(image_directory, exist_ok=True)

annotations_directory = os.path.join(base_directory, "annotations")
os.makedirs(annotations_directory, exist_ok=True)

# Instantiate a client, using the data portal GraphQL API by default
client = Client()

# Find all tomograms related to a specific organism
runs = Run.find(client, query_filters=[Run.dataset.id == 10000])

for run in tqdm(runs, desc="Runs"):
    dataset_name = run.name
    # make a folder for the annotations
    tomo_annotations_directory = os.path.join(
        annotations_directory, dataset_name
    )
    for voxel_spacing in tqdm(
        run.tomogram_voxel_spacings, desc="Voxel Spacings"
    ):
        for tomo in tqdm(voxel_spacing.tomograms, desc="Tomograms"):
            # download the tomogram
            tomo.download_mrcfile(image_directory)

            # This downloads all segmentation masks and .json files,
            # but not the .ndjson files containing point annotations
            tomo.download_all_annotations(
                tomo_annotations_directory,
                format="mrc",
                shape="SegmentationMask",
            )

        # Download point annotations
        for annotation in tqdm(voxel_spacing.annotations, desc="Annotations"):
            # download the annotation
            for annotation_file in annotation.files:
                if annotation_file.shape_type == "OrientedPoint" or annotation_file.shape_type == "Point":
                    annotation_file.download(tomo_annotations_directory)