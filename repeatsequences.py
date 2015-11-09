from bitmap import BitMap

class RepeatSequences:
    """Class to aggregate both repeat sequences and repeat
    sequence bitmap into a single object."""

    def __init__(self, size, frame_min, frame_max):
        self.repeats = dict()
        self.bitmap = BitMap(size)
        self.frame_min = frame_min
        self.frame_max = frame_max
        return

    """Getters."""
    def get_repeat_dict(self):
        return self.repeats
    def get_frame_min(self):
        return self.frame_min
    def get_frame_max(self):
        return self.frame_max
    def get_bitmap(self):
        return self.bitmap
