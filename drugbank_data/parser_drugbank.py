import xml.sax

class StructureHandler(xml.sax.ContentHandler):
    def __init__(self, output_file, max_depth=3):
        self.depth = 0
        self.max_depth = max_depth
        self.path = []
        self.output_file = output_file

    def startElement(self, name, attrs):
        self.path.append(name)
        if self.depth < self.max_depth:
            self.output_file.write('  ' * self.depth + name + '\n')
        self.depth += 1

    def endElement(self, name):
        self.path.pop()
        self.depth -= 1

def parse_large_xml(xml_file_path, output_file_path, max_depth=3):
    with open(output_file_path, 'w') as output_file:
        handler = StructureHandler(output_file, max_depth)
        parser = xml.sax.make_parser()
        parser.setContentHandler(handler)
        parser.parse(xml_file_path)

# Usage
xml_file_path = '/Users/sanazkazeminia/Documents/LLM_Agent/drugbank_data/full database.xml'
output_file_path = '/Users/sanazkazeminia/Documents/LLM_Agent/drugbank_data/structure.txt'
parse_large_xml(xml_file_path, output_file_path)

print(f"XML structure has been saved to {output_file_path}")

