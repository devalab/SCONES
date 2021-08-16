import typing
import numpy as np

from utils.sample_key import SampleKey
from utils.get_aligned_position import get_aligned_position

class SymmetricClosure:
    def __init__(self):
        self.data = {}

    def augment(self, data: dict, invertfn) -> None:
        for sample in data.values():
            fwd_sample, rev_sample = invertfn(sample)
            fwd_id, rev_id = fwd_sample["id"], rev_sample["id"]
            assert(fwd_id not in self.data)
            assert(rev_id not in self.data)
            self.data[fwd_id] = fwd_sample
            self.data[rev_id] = rev_sample

    def export(self, csv_path: str, columns: typing.Optional[list]=None) -> None:
        import pandas as pd
        df = pd.DataFrame.from_dict(self.data).T
        if columns:
            for column in columns:
                assert(column in df.columns)
            df = df.reindex(columns, axis=1)
            assert(df.columns.size == len(columns))
        df.to_csv(csv_path, index=False)

