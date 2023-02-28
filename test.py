from svgutils.compose import *
import os
import svgutils.transform as sg


import svg_stack as ss

doc = ss.Document()
svg_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.svg')]

layout1 = ss.HBoxLayout()
for ff in svg_files:
    layout1.addSVG(ff,ss.AlignHCenter)

doc.setLayout(layout1)
doc.save(os.path.join(input_dir,'final.svg'))


def merge_svg(input_dir, output_path):
    # Create a list of SVG files in the input directory
    
    # Load each SVG file as a SVG object
    svg_objs = [sg.fromfile(f).getroot() for f in svg_files]

    # Combine all SVG objects into one using the SVGFigure method
    combined_fig = sg.SVGFigure()
    for obj in svg_objs:
        combined_fig.append(obj)

    # Move all objects to (0, 0) to ensure they are visible
    combined_fig.move(0, 0)

    # Save the combined SVG to a file
    combined_fig.save(output_path)


input_dir = r"C:\Users\aswen\Documents\Data\Project_CRC\QC_temp\QCfigures"
output_path = r"C:\Users\aswen\Documents\Data\Project_CRC\QC_temp\QCfigures\final.svg"
merge_svg(input_dir, output_path)
