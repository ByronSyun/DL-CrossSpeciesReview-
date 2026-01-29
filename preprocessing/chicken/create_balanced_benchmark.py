import random

random.seed(42)

POSITIVE_FILE = "/mnt/userdata4/splicing/ChickenGTEx/processed_data/chicken_positive_samples_HIGH_QUALITY.tsv"
NEGATIVE_FILE = "/mnt/userdata4/splicing/ChickenGTEx/processed_data/chicken_negative_pool_HIGH_QUALITY.tsv"
OUTPUT_FILE = "/mnt/userdata4/splicing/ChickenGTEx/processed_data/chickengtex_silver_benchmark_balanced.tsv"
N_POSITIVE = 12500
N_NEGATIVE = 12500


def sample_file(input_file, n_samples, label):
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    if n_samples > len(lines):
        n_samples = len(lines)
    
    sampled = random.sample(lines, n_samples)
    processed = []
    for line in sampled:
        parts = line.strip().split('\t')
        if len(parts) >= 3:
            processed.append(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{label}\n")
    return processed


def main():
    positive_samples = sample_file(POSITIVE_FILE, N_POSITIVE, label=1)
    negative_samples = sample_file(NEGATIVE_FILE, N_NEGATIVE, label=0)
    all_samples = positive_samples + negative_samples
    random.shuffle(all_samples)
    
    with open(OUTPUT_FILE, 'w') as f:
        f.writelines(all_samples)
    
    print(f"Created {OUTPUT_FILE} with {len(all_samples)} samples")


if __name__ == "__main__":
    main()
