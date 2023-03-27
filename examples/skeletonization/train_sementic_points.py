import pytorch_lightning as pl
from monai.data import DataLoader
from monai.transforms import Compose, SpatialPadd
from pytorch_lightning.callbacks import ModelCheckpoint
from pytorch_lightning.loggers import TensorBoardLogger

from morphospaces.datasets import StandardHDF5Dataset
from morphospaces.networks.skeletonization import SemanticSkeletonNet
from morphospaces.transforms.image import (
    ExpandDimsd,
    ImageAsFloat32,
    MaskFromVectorField,
)

if __name__ == "__main__":
    batch_size = 1
    patch_shape = (96, 96, 96)
    patch_stride = (48, 48, 48)
    patch_threshold = 0.0005
    lr = 0.0002
    logdir_path = "./points_checkpoints"

    pl.seed_everything(42, workers=True)

    train_transform = Compose(
        [
            ImageAsFloat32(keys="raw"),
            MaskFromVectorField(input_key="raw", output_key="mask"),
            ExpandDimsd(keys="weights"),
            SpatialPadd(
                keys=["raw", "label", "weights", "mask"],
                spatial_size=[96, 96, 96],
            ),
        ]
    )

    train_ds = StandardHDF5Dataset.from_glob_pattern(
        glob_pattern="./test/*.h5",
        stage="train",
        transform=train_transform,
        patch_shape=patch_shape,
        stride_shape=patch_stride,
        patch_filter_ignore_index=(0,),
        patch_threshold=patch_threshold,
        patch_slack_acceptance=0.01,
        mirror_padding=(16, 32, 32),
        raw_internal_path="background_vector_image",
        label_internal_path="point_vectors",
        weight_internal_path="points_semantic_target",
    )
    train_loader = DataLoader(
        train_ds, batch_size=batch_size, shuffle=True, num_workers=4
    )

    val_transform = Compose(
        [
            ImageAsFloat32(keys="raw"),
            MaskFromVectorField(input_key="raw", output_key="mask"),
            ExpandDimsd(keys="weights"),
            SpatialPadd(
                keys=["raw", "label", "weights", "mask"],
                spatial_size=[96, 96, 96],
            ),
        ]
    )

    val_ds = StandardHDF5Dataset.from_glob_pattern(
        glob_pattern="./test/*.h5",
        stage="val",
        transform=val_transform,
        patch_shape=patch_shape,
        stride_shape=patch_stride,
        patch_filter_ignore_index=(0,),
        patch_threshold=patch_threshold,
        patch_slack_acceptance=0.005,
        mirror_padding=(16, 32, 32),
        raw_internal_path="background_vector_image",
        label_internal_path="point_vectors",
        weight_internal_path="points_semantic_target",
    )
    val_loader = DataLoader(
        val_ds, batch_size=batch_size, shuffle=False, num_workers=4
    )

    # make the checkpoint callback
    best_checkpoint_callback = ModelCheckpoint(
        save_top_k=1,
        monitor="val_loss",
        mode="min",
        dirpath=logdir_path,
        every_n_epochs=1,
        filename="vit-best",
    )
    last_checkpoint_callback = ModelCheckpoint(
        save_top_k=1,
        save_last=True,
        dirpath=logdir_path,
        every_n_epochs=1,
        filename="vit-last",
    )

    net = SemanticSkeletonNet(
        n_vectors_channels=6,
        out_channels=4,
        image_key="raw",
        vectors_gt_key="label",
        skeleton_gt_key="weights",
        learning_rate=lr,
    )

    # logger
    logger = TensorBoardLogger(
        save_dir=logdir_path, version=1, name="lightning_logs"
    )

    trainer = pl.Trainer(
        accelerator="gpu",
        devices=1,
        callbacks=[best_checkpoint_callback, last_checkpoint_callback],
        logger=logger,
        max_epochs=2000,
    )
    trainer.fit(
        net,
        train_dataloaders=val_loader,
        val_dataloaders=val_loader,
    )
