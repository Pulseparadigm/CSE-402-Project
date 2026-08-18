import pymupdf, glob

files = glob.glob('f:/4-1/402/project/*.pdf')
print("Found:", files)
doc = pymupdf.open(files[0])
for i, page in enumerate(doc):
    print(f"=== PAGE {i+1} ===")
    print(page.get_text())
