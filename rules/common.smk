import pandas as pd

# Helper function to import tables
def safely_read_csv(path, **kwargs):
    try:
        return pd.read_csv(path, **kwargs)
    except pd.errors.EmptyDataError:
        pass


# Helper function to concatenate output tables
def concatenate_tables(input, sep="\s+", cols_to_integer=None):
    frames = [safely_read_csv(f, sep=sep) for f in input]
    frames_concatenated = pd.concat(frames, keys=input, sort=False)
    if cols_to_integer:
        frames_concatenated[cols_to_integer] = frames_concatenated[
            cols_to_integer
        ].apply(lambda x: pd.Series(x, dtype="Int64"))
    return frames_concatenated
