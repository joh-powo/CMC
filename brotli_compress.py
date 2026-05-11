import brotli
import sys

def compress_brotli(input_file, output_file):
    with open(input_file, 'rb') as f:
        data = f.read()
    
    # Use maximum compression (quality=11)
    compressed = brotli.compress(data, quality=11)
    
    with open(output_file, 'wb') as f:
        f.write(compressed)
    
    print(f'Original: {len(data)} bytes')
    print(f'Compressed: {len(compressed)} bytes')
    print(f'Ratio: {(1 - len(compressed)/len(data))*100:.2f}%')

if __name__ == '__main__':
    if len(sys.argv) >= 3:
        compress_brotli(sys.argv[1], sys.argv[2])
    else:
        print('Usage: python brotli_compress.py input output')
