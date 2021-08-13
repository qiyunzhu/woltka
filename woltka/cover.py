class SortedRangeList:
    def __init__(self, autocompress=100000):
        self.ranges = []
        self.autocompress = autocompress
        self.since_compressed = 0

    def add_range(self, start, end):
        self.ranges.append((start, end))
        self.since_compressed += 1
        if self.autocompress is not None and \
                self.since_compressed > self.autocompress:
            self.compress()

    def compress(self):
        if self.since_compressed == 0:
            return
        # Sort ranges by start index
        self.ranges.sort(key=lambda r: r[0])

        new_ranges = []
        start_val = None
        end_val = None

        for r in self.ranges:
            if end_val is None:
                # case 1: no active range, start active range.
                start_val = r[0]
                end_val = r[1]
            elif end_val >= r[0] - 1:
                # case 2: active range continues through this range
                # extend active range
                end_val = max(end_val, r[1])
            else:  # if end_val < r[0] - 1:
                # case 3: active range ends before this range begins
                # write new range out, then start new active range
                new_range = (start_val, end_val)
                new_ranges.append(new_range)
                start_val = r[0]
                end_val = r[1]

        if end_val is not None:
            new_range = (start_val, end_val)
            new_ranges.append(new_range)

        self.ranges = new_ranges
        self.since_compressed = 0

    def compute_length(self):
        if self.since_compressed > 0:
            self.compress()
        total = 0
        for r in self.ranges:
            total += r[1] - r[0] + 1
        return total

    def __str__(self):
        return str(self.ranges)
