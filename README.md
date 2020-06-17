# BNDM-Analysis

An implementation of the classic BNDM algorithm,
modified based on Navarro and Raffinot's
[Fast and Flexible String Matching](https://dl.acm.org/doi/10.1145/351827.384246)
paper to support k-approximate matching. A wrapper
function (align) is used to implement a binary
search to find the least-error match and provide
an optimal alignment. Smith-Waterman and
Needleman-Wunsch are also implemented (based on
[this](https://tiefenauer.github.io/blog/smith-waterman/)
blog post and [this](https://gist.github.com/radaniba/11019717)
github repo) for comparison. A simple trial for
comparing the runtime of the three algorithms
on random text and a bank of dna sequences is
also implemented.

## bndm.py

uline() and hlight() simply wrap a string to
underline and highlight it respectively.

reconstruct() is a helper function used by bndm()
to generate an aligment from a match.

bndm() finds a match and generates the alignment.
@returns: remainder of text starting from match

align() finds the optimal alignment of the entire
pattern to a portion of the text with as few
errors as possible (with a log fraction threshold
on the number of allowed errors).
@returns: time taken in ns

Verbosity settings can be set to 0) no output,
1) output result, 2) display visualization.

## sw.py

matrix() and traceback() are internal functions
used for sw() to compute the cost matrix and find
the resulting optimal alignment respectively.

sw() is a wrapper function that gets the optimal
alignment of any portion of each text using the
Smith-Waterman algorithm.
@returns: time taken in ns, number of errors

## nw.py

matrix() and traceback() are internal functions
used for nw() to compute the cost matrix and find
the resulting optimal alignment respectively.

nw() is a wrapper function that gets the optimal
alignment of the entirety of each text using the
Needleman-Wunsch algorithm.
@returns: time taken in ns, number of errors

## analysis.py

generate() generates a random string from a
provided alphabet.

rand_letter() is a wrapper of generate() that
generates a single random character from a
provided alphabet.

get_similar() selects a random substring of the
text and generates random errors with a given
error probability (substitutions, insertions,
and deletions).

analyze() compares the runtime of a chosen
algorithm against the length of the pattern
searched at various text sizes.

load() converts the DNA sequences stored in
dna.txt into an array and strips metadata.

rand_seq() selects a random substring from a
random sequence in an array of sequences.

search_dna() compares the runtime of a chosen
algorithm against the length of the pattern
searched over various DNA sequences.

## Results

align() finds an optimal local-global alignment
between a pattern and a text, i.e. it aligns
the entire pattern to a portion of the text.
This is somewhere between the global
functionality of Needleman-Wunsch and the local
functionality of Smith-Waterman, both of which
require O(nm) space and time. align() gets away
with only O(mk) space where k is the maximum
number of errors allowed, which is a
significant improvement when the number of
errors is appropriately limited. The time
complexity is much more complex. When the
number of errors is sufficiently limited, the
runtime is approximately O(nlogm) since bndm()
takes approximately O(nklogm/m) time. However,
when enough errors are allowed, bndm() takes
closer to O(nm) time and the runtime of align()
approaches O(nm^2).

## Visualizations

To see the BNDM algorithm in action, check out
[this](https://bit.ly/2N2jmDP) link or the
included demo.gif file! This visualization was
rendered from demo.yml and was created using
[Terminalizer](https://terminalizer.com/)
with the included config.yml file.