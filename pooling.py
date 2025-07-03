from shiny import App, ui, reactive, render
import pandas as pd
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt

# Define the UI
app_ui = ui.page_fluid(
    #Centered Title
    ui.h2("Pooling", style="text-align: center;"),
    
    # Row with two input fields side by side (one on the left, one on the right)
    ui.row(
        ui.column(3),
        ui.column(
            5,
            ui.input_numeric("n_samp", "Number of Samples:", value=1)
        ),
        ui.column(
            2,
            ui.input_numeric("differentiate", "Max number of positives:", value=1)
        )
    ),
    
    # Submit button centered below inputs
    ui.row(
        ui.column(
            12,
            ui.input_action_button("submit", "Enter", style="display: block; margin-left: auto; margin-right: auto;")
        )
    ),
    
    # Text output for displaying last submitted values, centered below inputs
    ui.hr(),
    ui.div(
        ui.h4("Last submitted values:"),
        ui.output_text_verbatim("last_val"),
        style="text-align: center;"
    ),

    ui.hr(),
    ui.div(
        ui.h4("Database reply"),
        ui.output_text_verbatim("database_r"),  
        ui.div(
        ui.download_button("download_pickle", "Download complete pickled file"),
        style="text-align: center; margin-top: 60px;"
        ),
        style="text-align: center;"
    ),

    ui.hr(),
    ui.div(
        ui.h4("Summary pooling strategy table", style="text-align: center;"),
        # Wrap the table in a div and apply margin-left and margin-right for centering
        ui.div(
            ui.output_data_frame("summary_t"),
            style="display: flex; justify-content: center;"  # Flexbox centering
        ),
        ui.div(
        ui.download_button("download_summary", "Download summary"),
        style="text-align: center; margin-top: 60px;"
        )
    ),

    ui.hr(),
    ui.div(
        ui.h4("Command to run code locally"),
        ui.output_text_verbatim("commands"),  
        
        style="text-align: center;"
    ),


    ui.hr(),
    ui.div(
        ui.h4("Downloadable tables"),
        ui.output_text_verbatim("table_text"),  
        ui.download_button("download_table_matrix", "Download matrix pooling Table"),
        ui.download_button("download_table_md", "Download multidimensional pooling Table"),
        ui.download_button("download_table_random", "Download random pooling Table"),
        ui.download_button("download_table_STD", "Download STD pooling Table"),
        ui.download_button("download_table_CT", "Download Chinese trick pooling Table"),
        ui.download_button("download_table_binary", "Download binary pooling Table"),
        style="text-align: center;",
    ),


    # Conditionally display the summary table if extra_computation is False
    ui.hr(),
    ui.panel_conditional(
        "output.extra_computation === false",  # JavaScript condition to check if extra_computation is False
        ui.h4("Histogram Plot"),
        ui.output_plot("histogram_plot")  # Output container for the plot
    ),
    
    # Section for downloadable tables with dynamic names
    ui.hr(),
    ui.panel_conditional(
        "output.extra_computation === false",  # JavaScript condition to check if extra_computation is False
        ui.h4("Downloadable Tables"),
        ui.div(
            {"id": "download-section"},
            ui.output_ui("download_buttons")  # Placeholder for dynamically generated download buttons with names
        )
    ),
    
    # JavaScript for copying text to clipboard
    ui.tags.script("""
        function copyCommand(command) {
            navigator.clipboard.writeText(command).then(function() {
                alert('Copied to clipboard: ' + command);
            }).catch(function(error) {
                console.error('Failed to copy text: ', error);
            });
        }
    """),

    # JavaScript to trigger button click on Enter key press for both input fields
    ui.tags.script('''
        document.getElementById('n_samp').addEventListener('keydown', function(event) {
            if (event.key === 'Enter') {
                event.preventDefault();  // Prevent default form submission
                document.getElementById('submit').click();  // Trigger button click
            }
        });

        document.getElementById('differentiate').addEventListener('keydown', function(event) {
            if (event.key === 'Enter') {
                event.preventDefault();  // Prevent default form submission
                document.getElementById('submit').click();  // Trigger button click
            }
        });
    ''')
)
WA_DIRECTORY='D:\precomputed'
SCRAMBLER_DIRECTORY='D:\output'
MAX_DIFFERENTIATE=4

def find_n_folder(n_samp, wa_directory):
    folders = [f for f in os.listdir(wa_directory) if os.path.isdir(os.path.join(wa_directory, f)) and f.startswith('N_')]
    x_values = [int(f.split('_')[1]) for f in folders]
    if n_samp in x_values:
        return f'N_{n_samp}'
    greater_x = [x for x in x_values if x > n_samp]
    if greater_x:
        smallest_x = min(greater_x)
        return f'N_{smallest_x}'
    return None

def find_closest_diff_folder(n_folder_path, differentiate):
    folders = [f for f in os.listdir(n_folder_path) if os.path.isdir(os.path.join(n_folder_path, f)) and f.startswith('diff_')]
    y_values = [int(f.split('_')[1]) for f in folders]
    if differentiate in y_values:
        return f'diff_{differentiate}'
    if y_values:
        closest_y = min(y_values, key=lambda y: abs(y - differentiate))
        return f'diff_{closest_y}'
    return None

import os
import pandas as pd

def load_wa_matrices(folder_path):
    DFFS = {}
    # List all files in the folder
    files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]
    
    for file in files:
        # Check if file matches the WA_Method_N_x_diff_y pattern
        if file.startswith('WA_') and file.endswith('.xlsx'):
            # Example filename: WA_Method_N_10_diff_2.xlsx
            parts = file.split('_')
            # Extract method name (between WA_ and N_x)
            try:
                n_index = next(i for i, part in enumerate(parts) if part.startswith('N'))
                method_name = '_'.join(parts[1:n_index])
                # Load the Excel file as a DataFrame
                file_path = os.path.join(folder_path, file)
                matrix_df = pd.read_excel(file_path, header=None)
                # Rename columns and index
                matrix_df.columns = ['Pool ' + str(i) for i in range(matrix_df.shape[1])]
                matrix_df.index = ['Sample ' + str(i) for i in range(matrix_df.shape[0])]
                # Store in dictionary with method name as key
                DFFS[method_name] = matrix_df
            except StopIteration:
                # If no N_x part found, skip this file
                continue
    return DFFS




# Define the server logic
def server(input, output, session):

    output.last_values = reactive.Value("")
    output.database_reply = reactive.Value("")
    output.extra_computation = reactive.Value(0)
    output.personalized_command = reactive.Value("")
    output.dataframes = reactive.Value(0)
    output.full_pickle = reactive.Value(0)


    ls_met=['Pooling strategy', 'mean_experiments', 'max_compounds_per_well', 'n_wells', 'percentage_check', 'mean_extra_exp']
    output.summary_table = reactive.Value(pd.DataFrame(columns=ls_met))


    @reactive.Effect
    @reactive.event(input.submit)
    def _():
        # Get user inputs
        n_samp = input.n_samp()
        differentiate = input.differentiate()
        #print(n_samp)

        last_values_text = f"Max number of samples: {n_samp}, Max positives: {differentiate}"
        
        # Set output to display last submitted values
        output.last_values.set(last_values_text)

        if differentiate > n_samp:
            output_text=f'Maximum number of positives ({differentiate}) must always be smaller than the total number of samples ({n_samp})'
            output.database_reply.set(output_text)
            output.extra_computation.set(1)

        elif differentiate > MAX_DIFFERENTIATE:
            output_text=f'Maximum number of positives ({differentiate}) too high. The precomputed maximum is 4. To locally run the code for your specific setting follow the section below'
            output.database_reply.set(output_text)
            output.extra_computation.set(1)

        else:
            n_folder = find_n_folder(n_samp, WA_DIRECTORY)
            if n_folder:
                n_folder_path = os.path.join(WA_DIRECTORY, n_folder)
                diff_folder = find_closest_diff_folder(n_folder_path, differentiate)
                if diff_folder:
                    diff_folder_path = os.path.join(n_folder_path, diff_folder)
                    excel_filename = f'Metrics_{n_folder}_diff_{diff_folder.split("_")[1]}.xlsx'
                    excel_path = os.path.join(diff_folder_path, excel_filename)
                    if os.path.isfile(excel_path):
                        metrics_data = pd.read_excel(excel_path)

                        # Use metrics_data as needed
                    else:
                        # Handle missing metrics file
                        pass
                else:
                    # Handle missing diff_y folders
                    pass
            else:
                # Handle missing N_x folders
                pass

            '''
            full_dir='Final_precomputed_file.pk'
            with open(full_dir, 'rb') as handle:
                f1=pickle.load(handle)
            f2=f1['Differentiate '+str(differentiate)]
            a2=np.array(list(f2))
            md=np.max(a2)
            if n_samp>md:
                output_text=f'Maximum number of samples ({n_samp}) too high for the chosen number of positives \n Below are displayed the information for {md} samples and up to {differentiate} positives. \n You can run the code following the commands below for your specific case.'
                output.database_reply.set(output_text)
                output.extra_computation.set(0)

            elif np.sum(a2==n_samp)==0:
                md=np.min(a2[a2>n_samp])
                output_text=f'There is no precomputed strategy for {n_samp} samples. \n The closest precomputed strategy is for {md} samples with up to {differentiate} positives'
                output.database_reply.set(output_text)
                output.extra_computation.set(0)

            else:
                md=n_samp
                output_text=f'There is a precomputed strategy for {n_samp} samples with up to {differentiate} positives'
                output.database_reply.set(output_text)
                output.extra_computation.set(0)

            
            
            CR=f2[md]
            output.full_pickle.set(CR)
            '''

            #DFT=CR[0]
            #DFT.insert(loc=0, column='Pooling strategy', value=DFT.index)
            output.summary_table.set(metrics_data)

            #TBLS=CR[1]
            #DFFS={}
            for idx in TBLS.keys():
                b1=TBLS[idx]
                tmp1=pd.DataFrame(b1, columns=['Pool '+ str(i) for i in range(b1.shape[1])], index=['Sample '+ str(i) for i in range(b1.shape[0])])
                DFFS.update({idx:tmp1})

            output.dataframes.set(DFFS)

            
            
            
        command_p=f'python pooling_comparison.py --start {n_samp} --stop {n_samp+1} --step {1} --differentiate {differentiate} --rand_guesses {50}'

        output.personalized_command.set(command_p)

            #md=np.max(np.array(list(f2)))



        #md=np.max(np.array(list(f2)))



        if False:

            if differentiate > md:
                output.database_reply.set('output_text')

            md=np.max(np.array(list(f1)))
            if differentiate > np.max(np.array(list(f1))):
                output_text=f'Number of samples too high for chosen max number of positives ({differentiate}) <br> max number of samples pre-computed for pooling {md}. <br>'
                output.database_reply.set(output_text)
                extra_computation=True

            full_dir='Final_precomputed_file.pk'
            with open(full_dir, 'rb') as handle:
                f1=pickle.load(handle)
            
            md=np.max(np.array(list(f1)))
            if differentiate > np.max(np.array(list(f1))):
                output_text=f'Differentiate too high, max diffrentiate pre-computed {md}. <br>'
                extra_computation=True


            else:
                nm='Differentiate '+ differentiate
                f2=f1[nm]
                l2=list(f2)
                a2=np.array(l2)
                if n_samp in l2:
                    n_samp_new=n_samp.copy()

                elif n_samp<np.max(a2):
                    n_samp_new=a2[a2>n_samp][0] 

                else:
                    extra_computation=True

                if n_samp_new in l2:
                    f3=f2[n_samp]
                    summary=f3[0]
                    table_names=summary.index
                    tables=[i for i in f3.values()]



            download_buttons_html = "".join([
                f"<div><strong>{table_names[i]}</strong>: "
                f"<button id='download_table_{i}' onclick='Shiny.download(\"download_table_{i}\")'>Download {table_names[i]}</button></div><br>"
                for i in range(len(tables))
            ])

            output.download_buttons.set(ui.HTML(download_buttons_html))




            # Prepare terminal commands based on extra_computation results and user inputs
            commands = [
                f"echo 'git clone https://github.com/trouillon-lab/pooling.git'",  
                f"echo 'conda env create -n pooling --file=environments.yml'",
                f"echo 'conda activate pooling'",
                f"echo 'python pre-computation.py --start {n_samp} --stop {n_samp+1} --step 1 --differentiate {differentiate} --rand_guesses 10'" 
                f"process_data --samples {n_samp} --diff {differentiate}"
            ]

            # Generate HTML content for displaying commands with copy buttons
            command_html = "".join([
                f"<div><code>{command}</code> "
                f"<button onclick=\"copyCommand('{command}')\">Copy Command {i+1}</button></div><br>"
                for i, command in enumerate(commands)
            ])

            # Set the dynamically generated HTML content in the output UI
            output.commands_output.set(ui.HTML(command_html))


            # Set extra_computation flag to control conditional panel for plot
            output.extra_computation.set(extra_computation)


    @output
    @render.plot
    def histogram_plot():
        # Check if extra_computation is False before generating the plot
        #if not output.extra_computation.get():
        # Generate random data for histogram
        np.random.seed(19680801)
        data = 100 + 15 * np.random.randn(437)
        
        # Create a histogram with a fixed number of bins (or customize based on inputs)
        fig, ax = plt.subplots()
        ax.hist(data, bins=30, density=True)  # Fixed bins for simplicity or customize as needed
        ax.set_title("Histogram of Random Data")
        ax.set_xlabel("Value")
        ax.set_ylabel("Density")
        
        return fig

    @output
    @render.text
    def last_val():
        return  output.last_values.get()
    
    @output
    @render.text
    def database_r():
        return  output.database_reply.get()
    
    @output
    @render.text
    def commands():
        return  output.personalized_command.get()
    
    @output
    @render.data_frame
    def summary_t():      
        return output.summary_table.get()
    
    @output
    @render.download(filename=lambda: "summary.csv")
    async def download_summary():
        # Yield the content of the CSV file
        yield output.summary_table.get().to_csv(index=False)
    
    @output
    @render.download(filename=lambda: "matrix_pooling.csv")
    async def download_table_matrix():
        # Yield the content of the CSV file
        DFFS=output.dataframes.get()
        yield DFFS['matrix'].to_csv(index=True)
    
    @output
    @render.download(filename=lambda: "multidimensional_pooling.csv")
    async def download_table_md():
        # Yield the content of the CSV file
        DFFS=output.dataframes.get()
        idt=[i for i in list(DFFS) if i.startswith('multidim')]
        yield DFFS[idt[0]].to_csv(index=True)

    @output
    @render.download(filename=lambda: "random_pooling.csv")
    async def download_table_random():
        # Yield the content of the CSV file
        DFFS=output.dataframes.get()
        yield DFFS['random'].to_csv(index=True)
    
    @output
    @render.download(filename=lambda: "STD_pooling.csv")
    async def download_table_STD():
        # Yield the content of the CSV file
        DFFS=output.dataframes.get()
        yield DFFS['STD'].to_csv(index=True)

    @output
    @render.download(filename=lambda: "Chinese_trick_pooling.csv")
    async def download_table_CT():
        # Yield the content of the CSV file
        DFFS=output.dataframes.get()
        yield DFFS['Chinese trick'].to_csv(index=True)

    @output
    @render.download(filename=lambda: "binary_pooling.csv")
    async def download_table_binary():
        # Yield the content of the CSV file
        DFFS=output.dataframes.get()
        try:
            yield DFFS['Binary'].to_csv(index=True)
        except:
            yield pd.DataFrame().to_csv(index=True)

    @output
    @render.download(filename=lambda: "full_pooling.pk")
    async def download_pickle():
        with open('temp.pkl', 'wb') as f:
            pickle.dump(output.full_pickle, f)
            f.seek(0)  # Move back to start of file before yielding content
        
        # Read back and yield content for download
        with open('temp.pkl', 'rb') as f:
            yield f.read()


    #@output.plot("plot_placeholder")
    #def plot():
        # Placeholder for future plot implementation
    #    pass


# Create the app object
app = App(app_ui, server)