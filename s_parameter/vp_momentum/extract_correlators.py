#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
from compression import zstd
import json


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_filename")
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    parser.add_argument("--sink", action="append", dest="sinks")
    parser.add_argument("--source", action="append", dest="sources")
    parser.add_argument("--use_complex", action="store_true")
    parser.add_argument("--Nt", type=int, default=None)
    parser.add_argument("--Nx", type=int, default=None)
    parser.add_argument("--Ny", type=int, default=None)
    parser.add_argument("--Nz", type=int, default=None)
    return parser.parse_args()


def read_data(file_object, target_sinks, target_sources, use_complex=False):
    if use_complex:
        raise NotImplementedError("complex number support not yet implemented")

    data = {}
    mass = None

    reading_header_block = False
    reading_data_block = False

    for line in file_object:
        if line.startswith("start trajectory "):
            trajectory_index = int(line.split()[2])
            data[trajectory_index] = []
        if line.startswith("SOURCE:"):
            source = line.split()[1]
            if (not target_sources) or (source in target_sources):
                reading_header_block = True

        if line.startswith("MASSES:"):
            current_mass = float(line.split()[1])
            if mass is not None:
                if current_mass != mass:
                    raise ValueError("Inconsistent masses!")
            mass = current_mass

        if reading_header_block:
            if line.startswith("END"):
                reading_header_block = False
                reading_data_block = False
            if line.startswith("SOURCE_POS:"):
                source_position = tuple(map(int, line.split()[1:]))
                if (
                    not data
                    or not data[trajectory_index]
                    or data[trajectory_index][-1]["_source_position"] != source_position
                    or data[trajectory_index][-1]["_source"] != source
                ):
                    data[trajectory_index].append(
                        {"_source_position": source_position, "_source": source}
                    )

            if line.startswith("SINKS:") or line.startswith("SINK:"):
                sinks = line.split()[1:]
                if any(sink in sinks for sink in target_sinks):
                    reading_data_block = True
                    reading_header_block = False
                    slice_count = 0
                    datum = {sink: [] for sink in sinks}
                    continue

        if reading_data_block:
            if line.startswith("END"):
                reading_header_block = False
                reading_data_block = False
                data[trajectory_index][-1].update(datum)
                continue
            if line.startswith("SINK:"):
                data[trajectory_index][-1].update(datum)
                sink = line.split()[1]
                sinks = [sink]
                if sink not in target_sinks:
                    reading_data_block = False
                    reading_header_block = True
                else:
                    slice_count = 0
                    datum = {sink: []}
                continue

            split_line = line.split()
            assert int(split_line[0]) == slice_count

            if sinks == ["1LPBP"] and slice_count > 0:
                raise ValueError(f"t={slice_count} not expected for 1LPBP channel")

            slice_count += 1

            for sink_index, sink in enumerate(sinks):
                datum[sink].append(float(split_line[sink_index * 2 + 1]))

    return {
        "data": data,
        "mass": mass,
    }


def main():
    args = get_args()
    with zstd.open(args.input_filename, "rt") as file_object:
        data = read_data(file_object, args.sinks, args.sources, args.use_complex)

    # TODO: Check consistency between correlator length and given Nt
    data.update({"Nt": args.Nt, "Nx": args.Nx, "Ny": args.Ny, "Nz": args.Nz})
    json.dump(data, args.output_file)


if __name__ == "__main__":
    main()
