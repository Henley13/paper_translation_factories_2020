import config


class NucleiConfig(config.Config):
    NAME = "nuclei"
    GPU_COUNT = 1
    IMAGES_PER_GPU = 1
    NUM_CLASSES = 1 + 1 # background + nucleus
    TRAIN_ROIS_PER_IMAGE = 512
    STEPS_PER_EPOCH = 5000 # check mask_train for the final value
    VALIDATION_STEPS = 50
    DETECTION_MAX_INSTANCES = 512
    DETECTION_MIN_CONFIDENCE = 0.5
    DETECTION_NMS_THRESHOLD = 0.35
    RPN_NMS_THRESHOLD = 0.55
