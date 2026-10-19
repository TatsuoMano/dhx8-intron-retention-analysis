import sys
import csv

def parse_attributes(attr_str):
    attributes = {}
    for attr in attr_str.strip().split(';'):
        if not attr.strip():
            continue
        parts = attr.strip().split(' ')
        key = parts[0]
        value = parts[1].strip('"')
        attributes[key] = value
    return attributes

def format_attributes(attributes):
    return '; '.join(['{} "{}"'.format(k, v) for k, v in attributes.items()]) + ';'

def generate_introns(input_gtf, output_gtf):
    transcripts = {}
    
    print("Reading {}...".format(input_gtf))
    with open(input_gtf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            feature_type = parts[2]
            if feature_type != 'exon':
                continue
                
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parse_attributes(parts[8])
            
            if 'transcript_id' not in attributes:
                continue
            
            transcript_id = attributes['transcript_id']
            
            if transcript_id not in transcripts:
                transcripts[transcript_id] = {
                    'chrom': chrom,
                    'strand': strand,
                    'exons': [],
                    'attributes': attributes
                }
            
            transcripts[transcript_id]['exons'].append((start, end))

    print("Found {} transcripts. Generating introns...".format(len(transcripts)))
    
    with open(output_gtf, 'w') as out:
        for transcript_id, data in transcripts.items():
            exons = sorted(data['exons'], key=lambda x: x[0])
            
            if len(exons) < 2:
                continue
                
            for i in range(len(exons) - 1):
                intron_start = exons[i][1] + 1
                intron_end = exons[i+1][0] - 1
                
                if intron_start > intron_end:
                    continue
                
                # Create intron entry
                attr = data['attributes'].copy()
                # Remove exon specific attributes if any
                if 'exon_id' in attr: del attr['exon_id']
                if 'exon_number' in attr: del attr['exon_number']
                
                attr_str = format_attributes(attr)
                
                # GTF format: seqname source feature start end score strand frame attribute
                line = "{}\tgenerated\tintron\t{}\t{}\t.\t{}\t.\t{}\n".format(data['chrom'], intron_start, intron_end, data['strand'], attr_str)
                out.write(line)

    print("Done. Introns written to {}".format(output_gtf))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python generate_intron_gtf.py <input_gtf> <output_gtf>")
        sys.exit(1)
    
    generate_introns(sys.argv[1], sys.argv[2])
