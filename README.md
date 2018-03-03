# oxog-filter-legacy

Legacy filter for removing OxoG artifact (C>A/G>T in the CCG context), which was
originally created by Chip Stewart and Lee Lichtenstein.

It is assembled here for to serve as a benchmarking baseline.

## Use

### Option 1

Edit `oxog-filter.sh` to specify path to dependencies.
Then run `oxog-filter.sh` with required parameters.

### Option 2

Build docker image by

```
docker build -t oxog-filter-legacy .
```

Then the docker image may be run by

```
docker run --rm -v $(pwd):/root oxog-filter-legacy /bin/bash
```

Finally, run `oxog-filter.sh` with required parameters inside the docker shell.

