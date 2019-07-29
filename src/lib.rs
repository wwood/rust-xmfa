extern crate regex;

use std::io;
use std::fs;
use std::path::Path;
use std::io::BufRead;

use regex::Regex;
#[macro_use] extern crate lazy_static;

pub struct Reader<R: io::Read> {
    pub format_version: String,
    pub sequence_count: u64,
    pub sequence_file_names: Vec<String>,
    pub sequence_headers: Vec<String>,
    pub sequence_lengths: Vec<u64>,
    pub interval_count: u64,
    reader: io::BufReader<R>,
    line: String,
}

impl Reader<fs::File> {
    /// Read XMFA from given file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let f = fs::File::open(path)?;
        return Reader::new(f);
    }
}

impl <R: io::Read> Reader<R> {
    pub fn new(reader: R) -> Result<Self, io::Error> {
        let mut my_reader = io::BufReader::new(reader);

        let mut line_in = String::new();

        let mut format_version = None;
        let mut sequence_count = None;
        let mut sequence_file_names = vec!();
        let mut sequence_headers = vec!();
        let mut sequence_lengths = vec!();
        let mut interval_count = None;

        while my_reader.read_line(&mut line_in).is_ok() && line_in.get(0..1) == Some("#") {
            let line = line_in.trim_end().to_string();
            line_in = String::new();
            match line.find(" ") {
                None => panic!("Unexpected line in XMFA header: {}", line),
                Some(i) => {
                    match line.get(0..i) {
                        Some("#FormatVersion") => {
                            format_version = Some(
                                line.get((i+1)..)
                                    .expect(&format!(
                                        "Failed to parse FormatVersion header line: {}",
                                        line))
                                    .to_string());
                        },
                        Some("#SequenceCount") => {
                            sequence_count = Some(
                                line.get((i+1)..)
                                    .expect(&format!("Failed to parse SequenceCount header line: {}",
                                                     line))
                                    .parse::<u64>()
                                    .expect(&format!("Failed to convert SequenceCount line into an \
                                                      integer count: {}",
                                                     line)));
                        },
                        Some("##SequenceIndex") => {}, // Ignore for now.
                        Some("#IntervalCount") => {
                            interval_count = Some(
                                line.get((i+1)..)
                                    .expect(&format!("Failed to parse IntervalCount header line: {}",
                                                     line))
                                    .parse::<u64>()
                                    .expect(&format!("Failed to convert IntervalCount line into an \
                                                      integer count: {}",
                                                     line)));
                        },
                        Some("##SequenceFile") => {
                            sequence_file_names.push(
                                line.get((i+1)..)
                                    .expect(&format!("Failed to parse SequenceFile header line: {}",
                                                     line))
                                    .to_string()
                            )
                        },
                        Some("##SequenceHeader") => {
                            sequence_headers.push(
                                line.get((i+1)..)
                                    .expect(&format!("Failed to parse SequenceHeader header line: {}",
                                                     line)).to_string()
                            )
                        },
                        Some("##SequenceLength") => {
                            sequence_lengths.push(
                                line.get((i+1)..line.find("bp")
                                         .expect(&format!("Unexpectedly did not find 'bp' in SequenceLength line: {}", line)))
                                    .expect(&format!("Failed to parse SequenceLength header line: {}",
                                                     line))
                                    .parse::<u64>()
                                    .expect(&format!("Failed to convert SequenceLength line into an \
                                                      integer count: {}",
                                                     line))
                            )
                        },
                        _ => return Err(io::Error::new(
                            io::ErrorKind::Other, format!("Unexpected XMFA header entry: {}", line)))
                    }
                }
            }
        }

        Ok(Reader {
            format_version: format_version.expect("FormatVersion not found in header"),
            sequence_count: sequence_count.expect("SequenceCount not found in header"),
            sequence_file_names: sequence_file_names,
            sequence_headers: sequence_headers,
            sequence_lengths: sequence_lengths,
            interval_count: interval_count.expect("IntervalCount not found in header"),
            reader: my_reader,
            line: line_in,
        })
    }

    pub fn next_subalignment_block(&mut self) -> io::Result<(Vec<SubAlignment>)> {
        if self.line.get(0..1) != Some(">") {
            return Err(io::Error::new(
                io::ErrorKind::Other, format!(
                    "Unexpected XMFA sub-alignment block starting line: {}", self.line)))
        } else {
            let mut alignments = vec!();
            let mut current: Option<SubAlignment> = None;

            while self.line != "=\n" {
                //eprintln!("Line: {}", self.line);
                let line = self.line.trim_end();
                match line.get(0..1) {
                    Some(">") => {
                        if current.is_some() {
                            alignments.push(current.unwrap());
                        }

                        lazy_static! {
                            // >1:800-999 + cluster2 s1:p800
                            static ref RE: Regex = Regex::new(r"^>(\d+):(\d+)-(\d+) (.*)").unwrap();
                        }

                        match RE.captures(line) {
                            Some(mat) => {
                                current = Some(SubAlignment {
                                    sequence_number: mat.get(1).unwrap().as_str().parse::<u64>().unwrap(),
                                    start: mat.get(2).unwrap().as_str().parse::<u64>().unwrap(),
                                    stop: mat.get(3).unwrap().as_str().parse::<u64>().unwrap(),
                                    comment: mat.get(4).unwrap().as_str().to_string(),
                                    seq: String::new()
                                })
                            },
                            None => {
                                return Err(io::Error::new(
                                    io::ErrorKind::Other, format!(
                                        "Unexpected (type 2) XMFA sub-alignment block starting line: {}",
                                        self.line)))
                            }
                        }
                    },
                    _ => {
                        match current {
                            Some(ref mut a) => {a.seq.push_str(&line)},
                            _ => unreachable!()
                        }
                    }
                }
                self.line = String::new();
                self.reader.read_line(&mut self.line)?;
            }

            if current.is_some() {
                alignments.push(current.unwrap());
                self.line = String::new();
                self.reader.read_line(&mut self.line)?;
            }
            return Ok(alignments)
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct SubAlignment {
    pub sequence_number: u64,
    pub start: u64,
    pub stop: u64,
    pub comment: String,
    pub seq: String
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn it_reads_a_header() {
        let reader = Reader::from_file(Path::new("test/data/parsnp.xmfa")).unwrap();
        assert_eq!("Parsnp v1.1", reader.format_version);
        assert_eq!(3, reader.sequence_count);
        assert_eq!(
            vec![
                "1000_1.fna.ref",
                "1000_2.fna",
                "1000_3.fna",
            ], reader.sequence_file_names);
        assert_eq!(vec![
            ">genome1 random",
            ">genome2 one A to T SNP at position 100",
            ">genome3 G to C SNP at position 50, then take the last 200bp and put n front",
        ], reader.sequence_headers);
        assert_eq!(vec![1000,1000,1000], reader.sequence_lengths);
        assert_eq!(2, reader.interval_count);
    }

    #[test]
    fn it_reads_records() {
        let mut reader = Reader::from_file(Path::new("test/data/parsnp.xmfa")).unwrap();
        let mut aligns = reader.next_subalignment_block().unwrap();
        assert_eq!(vec!(
            SubAlignment {
                sequence_number: 1,
                start: 0,
                stop: 799,
                comment: "+ cluster1 :p1".to_string(),
                seq: "tcgcccgtcaccaccccaattcatacaccactagcggttagcaacgattGaccttgttatttatgccccgtgccctgaca\
                      gccttgaatggcttaccgtAatctatgagttagtaatttcgcctccttgatccaatcgtcgctcactccgcccggaaaac\
                      cagttaatcgacctggtgtatcgtacagaccacatctttctctgcataggcttagccggtgttccagaatcgttctcgag\
                      caatattcctggctgaatgaggaatacacgcgatctgtcaaatctggatctgcctacagtactaacgcaccagccggggt\
                      aactggggcgagccccctcggttacccgtggggtactacatgcttgagggtcgatactagccagcactgtactttcggtg\
                      ccaaggccactcgacgtatgttgtacgtgggacacgccgagatggactgggccaccaatattggaaggtttatacctgca\
                      tcatagagatacaaggaatactagccttgaacgtatggtaataattatcgctaggcccctatgctcttacaacccgttat\
                      tatgacccatacaactaataagtccgtttcattggttcagcagccctcatgtaccataccctcgttgcggatttatagcc\
                      agctctaggaccggggacattaccgcggaaaaacgcatttctttaggagagcgcttggacacccctgcgatcctcctatt\
                      ccttgcaagggacaacagttgcatagttcgcaactccattgacgggtatatacttttgaaacacgttggccgcgccacaa".to_string()
            },
            SubAlignment {
                sequence_number: 2,
                start: 0,
                stop: 799,
                comment: "+ cluster1 :p1".to_string(),
                seq: "tcgcccgtcaccaccccaattcatacaccactagcggttagcaacgattGaccttgttatttatgccccgtgccctgaca\
                      gccttgaatggcttaccgtTatctatgagttagtaatttcgcctccttgatccaatcgtcgctcactccgcccggaaaac\
                      cagttaatcgacctggtgtatcgtacagaccacatctttctctgcataggcttagccggtgttccagaatcgttctcgag\
                      caatattcctggctgaatgaggaatacacgcgatctgtcaaatctggatctgcctacagtactaacgcaccagccggggt\
                      aactggggcgagccccctcggttacccgtggggtactacatgcttgagggtcgatactagccagcactgtactttcggtg\
                      ccaaggccactcgacgtatgttgtacgtgggacacgccgagatggactgggccaccaatattggaaggtttatacctgca\
                      tcatagagatacaaggaatactagccttgaacgtatggtaataattatcgctaggcccctatgctcttacaacccgttat\
                      tatgacccatacaactaataagtccgtttcattggttcagcagccctcatgtaccataccctcgttgcggatttatagcc\
                      agctctaggaccggggacattaccgcggaaaaacgcatttctttaggagagcgcttggacacccctgcgatcctcctatt\
                      ccttgcaagggacaacagttgcatagttcgcaactccattgacgggtatatacttttgaaacacgttggccgcgccacaa".to_string()
            },
            SubAlignment {
                sequence_number: 3,
                start: 200,
                stop: 999,
                comment: "+ cluster1 s1:p200".to_string(),
                seq: "tcgcccgtcaccaccccaattcatacaccactagcggttagcaacgattCaccttgttatttatgccccgtgccctgaca\
                      gccttgaatggcttaccgtAatctatgagttagtaatttcgcctccttgatccaatcgtcgctcactccgcccggaaaac\
                      cagttaatcgacctggtgtatcgtacagaccacatctttctctgcataggcttagccggtgttccagaatcgttctcgag\
                      caatattcctggctgaatgaggaatacacgcgatctgtcaaatctggatctgcctacagtactaacgcaccagccggggt\
                      aactggggcgagccccctcggttacccgtggggtactacatgcttgagggtcgatactagccagcactgtactttcggtg\
                      ccaaggccactcgacgtatgttgtacgtgggacacgccgagatggactgggccaccaatattggaaggtttatacctgca\
                      tcatagagatacaaggaatactagccttgaacgtatggtaataattatcgctaggcccctatgctcttacaacccgttat\
                      tatgacccatacaactaataagtccgtttcattggttcagcagccctcatgtaccataccctcgttgcggatttatagcc\
                      agctctaggaccggggacattaccgcggaaaaacgcatttctttaggagagcgcttggacacccctgcgatcctcctatt\
                      ccttgcaagggacaacagttgcatagttcgcaactccattgacgggtatatacttttgaaacacgttggccgcgccacaa".to_string()
            }
        ), aligns);


        aligns = reader.next_subalignment_block().unwrap();
        assert_eq!(vec!(
            SubAlignment {
                sequence_number: 1,
                start: 800,
                stop: 999,
                comment: "+ cluster2 s1:p800".to_string(),
                seq: "atctcgactatcccgcatcagaacgcggtttggctccgaatgtacacaaataaccaggtatctatgatgtaatccggcga\
                      ctgtccagcttagaggcggtccttattacaccgaatggtaagcagatatgcgggccaacgaaaaatagttgtattgggcg\
                      cggccatggactttgcatcgagctgtccgtgtggcgaatc".to_string()
            },
            SubAlignment {
                sequence_number: 2,
                start: 800,
                stop: 999,
                comment: "+ cluster2 s1:p800".to_string(),
                seq: "atctcgactatcccgcatcagaacgcggtttggctccgaatgtacacaaataaccaggtatctatgatgtaatccggcga\
                      ctgtccagcttagaggcggtccttattacaccgaatggtaagcagatatgcgggccaacgaaaaatagttgtattgggcg\
                      cggccatggactttgcatcgagctgtccgtgtggcgaatc".to_string()
            },
            SubAlignment {
                sequence_number: 3,
                start: 0,
                stop: 199,
                comment: "+ cluster2 :p1".to_string(),
                seq: "atctcgactatcccgcatcagaacgcggtttggctccgaatgtacacaaataaccaggtatctatgatgtaatccggcga\
                      ctgtccagcttagaggcggtccttattacaccgaatggtaagcagatatgcgggccaacgaaaaatagttgtattgggcg\
                      cggccatggactttgcatcgagctgtccgtgtggcgaatc".to_string()
            },
        ), aligns);

        assert!(reader.next_subalignment_block().is_err());
    }
}
