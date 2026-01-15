
from shiny import App, ui, reactive, render
import pandas as pd
import os
import pickle
import numpy as np
import os
import pandas as pd
import scipy.stats
#import matplotlib.pyplot as plt
summary_text = ""
# Define the UI
app_ui = ui.page_fluid(
    ui.tags.title("PoolPy"),
    ui.navset_tab(
        ui.nav_panel("Pooling Methods Comparison",
            # Centered Title
            ui.h2("PoolPy", style="text-align: center;"),
            # Short description below the title
            ui.div(
                "This app evaluates combinatorial pooling performances of different methods. ",
                ui.br(),
                "For an overview head down to the ",
                ui.a("methods and metrics", href="#methods-metrics", id="nav-methods-link", style="text-decoration: underline; cursor: pointer;"),
                " section.",
                " For a detailed description read our ",
                ui.a("paper", href="#paper-link", style="text-decoration: underline;"),
                ".",
                ui.br(),
                "If you are working with sample prevalence and are unsure about the number of positives in your experiment, feel free to consult the ",
                ui.a("prevalence", href="#prevalence-section", id="nav-prevalence-link", style="text-decoration: underline; cursor: pointer;"),
                " section.",
                style="text-align: center; margin-bottom: 18px; font-size: 16px; color: #444;"
            ),
            ui.tags.script('''
                document.addEventListener('DOMContentLoaded', function() {
                    var methodsLink = document.getElementById('nav-methods-link');
                    if (methodsLink) {
                        methodsLink.addEventListener('click', function(event) {
                            event.preventDefault();
                            var tab = document.querySelector("[data-value='Methods & Metrics']");
                            if (tab) tab.click();
                        });
                    }
                    var prevalenceLink = document.getElementById('nav-prevalence-link');
                    if (prevalenceLink) {
                        prevalenceLink.addEventListener('click', function(event) {
                            event.preventDefault();
                            var tab = document.querySelector("[data-value='Prevalence']");
                            if (tab) tab.click();
                        });
                    }
                });
            '''),
            
            # Row with two input fields side by side (one on the left, one on the right)
            ui.row(
                ui.column(4),
                ui.column(2, ui.input_numeric("n_samp", "Number of Samples:", value=20)),
                ui.column(1),
                ui.column(2, ui.input_numeric("differentiate", "Max number of positives:", value=1)),
                ui.column(2),
            ),
            # ...existing code for Main page...
            ui.row(
                ui.column(
                    12,
                    ui.input_action_button("submit", "Enter", style="display: block; margin-left: auto; margin-right: auto;")
                )
            ),
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
                style="text-align: center;"
            ),
            ui.panel_conditional(
                "output.allow_downloads == 'true'",
                ui.div(
                    ui.h4("Summary pooling strategy table", style="text-align: center;"),
                    ui.div(
                        ui.output_ui("summary_colored_table"),
                        style="display: flex; justify-content: center;"
                    ),

# Add a new output to render the colored summary table in the server function, not in the UI definition
                    ui.div(
                        ui.output_ui("summary_text"),
                        style="margin: 18px 0 0 0; color: #555; font-size: 15px; text-align: center;"
                    ),
                    ui.div(
                        ui.download_button("download_summary", "Download summary"),
                        style="text-align: center; margin-top: 60px;"
                    )
                )
            ),
            ui.hr(),
            ui.div(
                ui.h4("Command to run code locally"),
                ui.output_text_verbatim("commands"),
                style="text-align: center;"
            ),
            ui.panel_conditional(
                "output.allow_downloads == 'true'",
                ui.div(
                    ui.h4("Downloadable tables"),
                    ui.download_button("download_table_matrix", "Matrix"),
                    ui.download_button("download_table_2d", "2-dim"),
                    ui.download_button("download_table_3d", "3-dim"),
                    ui.download_button("download_table_4d", "4-dim"),
                    ui.download_button("download_table_random", "Random"),
                    ui.download_button("download_table_STD", "STD"),
                    ui.download_button("download_table_CT", "Chinese remainder"),
                    ui.download_button("download_table_CTB", "Chinese rm bktrk"),
                    ui.download_button("download_table_CTS", "Chinese rm special"),
                    ui.download_button("download_table_binary", "Binary"),
                    style="text-align: center;",
                )
            ),
            ui.hr(),
            ui.panel_conditional(
                "output.allow_downloads == 'true'",
                ui.div(
                    ui.h4("Command to run decoder locally"),
                    ui.output_text_verbatim("decode"),
                    style="text-align: center;"
                ),
            ),
            ui.panel_conditional(
                "output.extra_computation == 'true'",
                ui.h4("Histogram Plot"),
                ui.output_plot("histogram_plot")
            ),
            ui.panel_conditional(
                "output.extra_computation == 'true'",
                ui.h4("Downloadable Tables"),
                ui.div(
                    {"id": "download-section"},
                    ui.output_ui("download_buttons")
                )
            ),
            ui.div(
                ui.output_text_verbatim("allow_downloads"),
                style="position: fixed; bottom: 2px; right: 2px; opacity: 0.05; font-size: 8px; z-index: 1; pointer-events: none;"
            ),
            ui.div("", style="height: 80px;"),
            ui.tags.script("""
                function copyCommand(command) {
                    navigator.clipboard.writeText(command).then(function() {
                        alert('Copied to clipboard: ' + command);
                    }).catch(function(error) {
                        console.error('Failed to copy text: ', error);
                    });
                }
            """),
            ui.tags.script('''
                document.getElementById('n_samp').addEventListener('keydown', function(event) {
                    if (event.key === 'Enter') {
                        event.preventDefault();
                        document.getElementById('submit').click();
                    }
                });
                document.getElementById('differentiate').addEventListener('keydown', function(event) {
                    if (event.key === 'Enter') {
                        event.preventDefault();
                        document.getElementById('submit').click();
                    }
                });
            ''')
        ),
        ui.nav_panel("Methods & Metrics",
            ui.h2("Methods & Metrics", style="text-align: center;"),
            ui.div(
                ui.h4("Pooling Methods Tested:"),
                ui.tags.ul(
                    ui.tags.li(ui.tags.b("Matrix"), ": Classic matrix pooling, each sample is included in a unique row and column pool."),
                    ui.tags.li(ui.tags.b("Multidimensional (2D, 3D, 4D)"), ": Samples are arranged in higher-dimensional grids, each dimension representing a pooling axis."),
                    ui.tags.li(ui.tags.b("Random"), ": Pools are formed by randomly assigning samples, often used as a baseline."),
                    ui.tags.li(ui.tags.b("STD"), ": Standard design, a reference pooling method with fixed structure."),
                    ui.tags.li(ui.tags.b("Chinese Remainder (bktrk, special)"), ": Pooling based on the Chinese Remainder Theorem, with variants for backtracking (bktrk) and special cases (special)."),
                    ui.tags.li(ui.tags.b("Binary"), ": Each sample is assigned to pools according to a binary code, maximizing information per test.")
                ),
                ui.h4("Metrics in the Summary Table:"),
                ui.tags.ul(
                    ui.tags.li(ui.tags.b("Pooling strategy"), ": Name of the pooling method used."),
                    ui.tags.li(ui.tags.b("Mean experiments"), ": Average number of tests required to identify the positive samples."),
                    ui.tags.li(ui.tags.b("Max compounds per well"), ": Maximum number of samples combined in any single pool."),
                    ui.tags.li(ui.tags.b("N pools"), ": Total number of pools used in first step of the strategy. (or in all steps for hierachical clustering)"),
                    ui.tags.li(ui.tags.b("Percentage check"), ": Fraction of cases requiring additional verification or retesting."),
                    ui.tags.li(ui.tags.b("Mean extra experiments"), ": Average number of extra tests needed beyond the initial pooling scheme.")
                ),
                style="max-width: 700px; margin: 0 auto; font-size: 16px; color: #333;"
            )
        ),
        ui.nav_panel("Prevalence",
                ui.h2("Prevalence", style="text-align: center;"),
                ui.div(
                    ui.div(
                        "This section is meant to guide the decision of the best pooling strategy for each individual case. Below you will find two numerical tables.",
                        ui.br(),
                        "The right table reports the probability of making at least one mistake while reading the results of one combinatorial pooling.",
                        ui.br(),
                        "The left table reports the probability of making at least one mistake while reading the results of all the combinatorial poolings of the experiment combined.",
                        style="text-align: center; color: #444; font-size: 16px;"
                    ),
                    style="text-align: center; margin-bottom: 18px;"
                ),
                ui.div(
                    ui.row(
                        ui.column(2),
                        ui.column(3,
                            ui.input_numeric("prev_n_samp", "Number of Samples:", value=100),
                        ),
                        ui.column(3,
                            ui.div(
                                ui.input_numeric(
                                    "prev_prevalence",
                                    "Prevalence (as fraction):",
                                    value=0.01,
                                ),
                                style="min-width: 220px; max-width: 100%; width: 100%; display: flex; align-items: center;"
                            )
                        ),
                        ui.column(3,
                            ui.input_numeric("prev_max_error", "Maximum acceptable error:", value=0.05),
                        ),
                    ),
                    style="margin-top: 30px;"
                ),
                ui.div(
                    ui.input_action_button("prev_submit", "Enter", style="display: block; margin: 20px auto 0 auto;"),
                    style="text-align: center;"
                ),
                ui.tags.script('''
                    document.addEventListener('DOMContentLoaded', function() {
                        ["prev_n_samp", "prev_prevalence", "prev_max_error"].forEach(function(id) {
                            var el = document.getElementById(id);
                            if (el) {
                                el.addEventListener('keydown', function(event) {
                                    if (event.key === 'Enter') {
                                        event.preventDefault();
                                        var btn = document.getElementById('prev_submit');
                                        if (btn) btn.click();
                                    }
                                });
                            }
                        });
                    });
                '''),
                ui.hr(),
                ui.div(
                    ui.h4("Last submitted values:"),
                    ui.output_text_verbatim("last_val_prev"),
                    style="text-align: center;"
                ),
                ui.div(
                    ui.output_ui("prevalence_error_table"),
                    style="margin-top: 30px; display: flex; justify-content: center;"
                ),
                ui.div(
                    ui.output_ui("prevalence_legend"),
                    style="margin-top: 18px; display: flex; justify-content: center;"
                ),
                ui.div(
                    ui.output_ui("prevalence_explanation"),
                    style="margin-top: 8px; display: flex; justify-content: center;"
                ),
            )
    )
)


WA_SUB_DIRECTORY='precomputed'
#SCRAMBLER_DIRECTORY='.\output'
MAX_DIFFERENTIATE=5
MAX_N=175

def find_n_folder(n_samp, wa_directory):
    folders = [f for f in os.listdir(wa_directory) if os.path.isdir(os.path.join(wa_directory, f)) and f.startswith('N_')]
    #print(folders)
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
    valid_folders = []
    valid_y_values = []
    for f in folders:
        y = int(f.split('_')[1])
        metrics_file = f"Metrics_{os.path.basename(n_folder_path)}_diff_{y}.csv"
        metrics_path = os.path.join(n_folder_path, f, metrics_file)
        if os.path.isfile(metrics_path):
            valid_folders.append(f)
            valid_y_values.append(y)
    if differentiate in valid_y_values:
        return f'diff_{differentiate}', differentiate
    elif valid_y_values:
        closest_y = min(valid_y_values, key=lambda y: abs(y - differentiate))
        return f'diff_{closest_y}', closest_y
    return None

def find_closest_n_folder(n, wa_directory):
    folders = [f for f in os.listdir(wa_directory) if os.path.isdir(os.path.join(wa_directory, f)) and f.startswith('N_')]
    x_values = [int(f.split('_')[1]) for f in folders]
    # Filter for x >= n
    valid_x = [x for x in x_values if x >= n]
    if valid_x:
        closest_x = min(valid_x)
        return f'N_{closest_x}', closest_x
    return None


# Function to generate expected error table for prevalence tab
def expected_error_table(N, prevalence, max_diff=None, correct=False):
    if max_diff is None:
        max_diff = 4  # or set to N or another reasonable value
    diff_values = np.arange(1, int(max_diff)+1)
    pop_sizes = []
    min_N=5
    N_current = int(np.ceil(N))
    factor=[]
    NF=1
    # Build the population sizes as described, always take the ceiling to be on the safe side
    while N_current >= min_N:
        pop_sizes.append(N_current)
        factor.append(NF)
        N_current = int(np.ceil(N_current / 2))
        NF=int(NF*2)
        if N_current < min_N:
            break
        pop_sizes.append(N_current)
        factor.append(NF)
        NF=int(NF*2.5)
        N_current = int(np.ceil(N_current / 2.5))
        if N_current < min_N:
            break
        pop_sizes.append(N_current)
        factor.append(NF)
        NF=int(NF*2)
        N_current = int(np.ceil(N_current / 2))
    # Remove duplicates and sort descending
    #pop_sizes = sorted(set([x for x in pop_sizes if x >= min_N]), reverse=True)
    
    # Build the table
    data = []
    for N_val,NF in zip(pop_sizes,factor):
        row = []
        for diff in diff_values:
            p_error = 1 - scipy.stats.binom.cdf(diff, N_val, prevalence)
            if correct:
                p_error=1-(1-p_error)**int(NF)
            row.append(p_error)
        data.append(row)
    df = pd.DataFrame(data, index=pop_sizes, columns=diff_values)
    # No minimum error logic needed
    df.index.name = "N"
    df.columns.name = "Differentiate"
    return df





def load_wa_matrices(folder_path):
    DFFS = {}
    # List all files in the folder
    files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]
    
    for file in files:
        # Check if file matches the WA_Method_N_x_diff_y pattern
        if file.startswith('WA_') and file.endswith('.csv'):
            # Example filename: WA_Method_N_10_diff_2.xlsx
            parts = file.split('_')
            # Extract method name (between WA_ and N_x)
            try:
                n_index = next(i for i, part in enumerate(parts) if part.startswith('N'))
                method_name = '_'.join(parts[1:n_index])
                # Load the Excel file as a DataFrame
                file_path = os.path.join(folder_path, file)
                matrix_df = pd.read_csv(file_path, header=None)
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

    @output
    @render.ui
    def summary_colored_table():
        df = output.summary_table.get()
        if df is None or df.empty:
            return ui.HTML("<div>No data available.</div>")
        # Only color numeric columns
        df_numeric = df.select_dtypes(include=[float, int])
        # Build HTML table
        html = '<table style="border-collapse: collapse; font-size: 16px;">'
        # Header
        html += '<tr>'
        for col in df.columns:
            th_style = 'border: 1px solid #aaa; padding: 8px 12px; text-align: center;'
            if col == 'N pools':
                th_style += ' min-width: 120px;'
            html += f'<th style="{th_style}">{col}</th>'
        html += '</tr>'

        # Color logic
        green_indices = {}
        # Green for three smallest in 'Mean experiments' and 'Max compounds in a pool'
        for col in ['Mean experiments', 'Max compounds in a pool']:
            if col in df.columns:
                col_vals = df[col]
                idxs = pd.Series(col_vals).nsmallest(3).index.tolist()
                green_indices[col] = idxs


        # Red or green for 'Mean steps' (if <=3 with >1, color those red; if >3, color all ==1 green)
        red_indices = []
        if 'Mean steps' in df.columns:
            col_vals = df['Mean steps']
            gt1 = col_vals[col_vals > 1]
            if 1 <= len(gt1) <= 3:
                red_indices = gt1.index.tolist()
            elif len(gt1) > 3:
                green_indices['Mean steps'] = col_vals.index[col_vals == 1].tolist()

        # Rows
        for i, row in df.iterrows():
            html += '<tr>'
            for col in df.columns:
                val = row[col]
                td_style = 'border: 1px solid #aaa; padding: 8px 12px; text-align: center;'
                if col == 'N pools':
                    td_style += ' min-width: 120px;'
                # Green logic
                if col in green_indices and i in green_indices[col]:
                    td_style += ' background-color: #b6fcb6; color: #111;'
                # Red logic for Mean steps
                if col == 'Mean steps' and i in red_indices:
                    td_style += ' background-color: #ffcccc; color: #111;'
                html += f'<td style="{td_style}">{val}</td>'
            html += '</tr>'
        html += '</table>'
        return ui.HTML(html)
    # For Prevalence tab: track last submitted values
    output.last_val_prev = reactive.Value("")
    # Prevalence error table output

    # Prevalence tab reactives
    output.prev_n_samp = reactive.Value(100)
    output.prev_prevalence = reactive.Value(0.01)
    output.prev_max_error = reactive.Value(0.05)

            
    output.summary_text = reactive.Value("")
    output.last_values = reactive.Value("")
    output.database_reply = reactive.Value("")
    output.extra_computation = reactive.Value(0)
    output.personalized_command = reactive.Value("")
    output.decoder = reactive.Value("")
    output.dataframes = reactive.Value(0)
    output.full_pickle = reactive.Value(0)
    output.debug = reactive.Value("")
    output.allow_downloads = reactive.Value(False)






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
            output.summary_table.set(pd.DataFrame(columns=ls_met))
            output.allow_downloads.set(False)

        elif differentiate > MAX_DIFFERENTIATE:
            output_text_1=f'Maximum number of positives ({differentiate}) too high. The precomputed maximum is {MAX_DIFFERENTIATE}.\n'
            output_text_2='To locally run the code for your specific setting follow the section below'
            output_text=output_text_1+output_text_2
            output.database_reply.set(output_text)
            output.extra_computation.set(1)
            output.summary_table.set(pd.DataFrame(columns=ls_met))
            output.allow_downloads.set(False)

        elif n_samp > MAX_N:
            output_text_1=f'Maximum number of samples ({n_samp}) too high. The precomputed maximum is {MAX_N}.\n'
            output_text_2='To locally run the code for your specific setting follow the section below'
            output_text=output_text_1+output_text_2
            output.database_reply.set(output_text)
            output.extra_computation.set(1)
            output.summary_table.set(pd.DataFrame(columns=ls_met))
            output.allow_downloads.set(False)

        else:
            try:
                app_root = os.listdir('.')[0]
                WA_DIRECTORY=os.path.join(app_root,WA_SUB_DIRECTORY)
                #n_folder = find_n_folder(n_samp, WA_DIRECTORY)
                n_folder, new_n = find_closest_n_folder(n_samp, WA_DIRECTORY)
                if n_folder:
                    n_folder_path = os.path.join(WA_DIRECTORY, n_folder)
                    diff_folder, new_diff = find_closest_diff_folder(n_folder_path, differentiate)
                    if diff_folder:
                        diff_folder_path = os.path.join(n_folder_path, diff_folder)
                        excel_filename = f'Metrics_{n_folder}_diff_{diff_folder.split("_")[1]}.csv'
                        excel_path = os.path.join(diff_folder_path, excel_filename)
                        output.allow_downloads.set(True)
                        if new_n!=n_samp or new_diff!=differentiate:
                            output_text_1=f'There is no precomputed strategy for {n_samp} samples with up to {differentiate} positives \n'
                            output_text_2=f'The next best precomputed strategy is for {new_n} samples with up to {new_diff} positives shown below'
                            output_text=output_text_1+output_text_2
                            output.database_reply.set(output_text)
                        else: 
                            output_text=f'There is a precomputed strategy for {n_samp} samples with up to {differentiate} positives'
                            output.database_reply.set(output_text)
                        #print(output_text)
                        if os.path.isfile(excel_path):
                            metrics_data = pd.read_csv(excel_path)

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
                    output.extra_computation.set(0)I go

                
                
                CR=f2[md]
                output.full_pickle.set(CR)
                '''

                #DFT=CR[0]
                #DFT.insert(loc=0, column='Pooling strategy', value=DFT.index)
                metrics_data.drop(metrics_data.columns[0], axis=1, inplace=True)
                #ls_met_nice=['Pooling design', 'Mean total experiments', 'Max samples per pool', 'N pools', '\% multiple rounds', 'Mean extra experiments']
                dict_ren={'N wells':'N pools', 'Max compunds per well': 'Max compounds in a pool'}#{i:j for i,j in zip(metrics_data.columns,ls_met_nice)}
                metrics_data.rename(columns=dict_ren, inplace=True)
                #metrics_data.rename(index={'Chinese trick':'Chinese reminder'}, inplace=True)
                output.summary_table.set(metrics_data)

                copio=metrics_data.copy()
                if 'Mean experiments' in copio.columns and 'Mean steps' in copio.columns:
                    copio = copio.sort_values(by=['Mean experiments', 'Mean steps'], ascending=[True, True])
                long_exps=f'If your experiments are <b>expensive</b>, the best method for you might be <b>{copio.iloc[0,0]}</b> one.<br>'
                
                if 'Mean experiments' in copio.columns and 'Mean steps' in copio.columns:
                    copio = copio.sort_values(by=['Mean steps', 'Mean experiments'], ascending=[True, True])
                
                expensive_exps=f'If your experiments are <b>long</b>, the best method for you might be the <b>{copio.iloc[0,0]}</b> one.<br>'

                if 'Mean experiments' in copio.columns and 'Max compounds in a pool' in copio.columns:
                    copio = copio.sort_values(by=['Max compounds in a pool', 'Mean experiments'], ascending=[True, True])                

                low_pools=f'If you <b>cannot pool many samples</b> together, the best method for you might be <b>{copio.iloc[0,0]}</b> one.<br>'

                output.summary_text.set(expensive_exps+long_exps+low_pools)

                #TBLS=CR[1]
                DFFS={}
                table_path=os.path.join(diff_folder_path,'WAs')
                TBLS=load_wa_matrices(table_path)
                #print(table_path)
                for idx in TBLS.keys():
                    b1=TBLS[idx]
                    tmp1=pd.DataFrame(b1, columns=['Pool '+ str(i) for i in range(b1.shape[1])], index=['Sample '+ str(i) for i in range(b1.shape[0])])
                    DFFS.update({idx:tmp1})

                output.dataframes.set(DFFS)
                #print(DFFS.keys())

            except Exception as e:
                #output.debug.set(str(e))
                files = os.listdir('.')
                app_root = os.listdir('.')[0]
                output.debug.set(str([f for f in files]))
                output.debug.set(f"cwd: {os.getcwd()}, files: {os.listdir('.')}")

                app_root = os.listdir('.')[0]  # first (and usually only) item — the app folder
                app_root_contents = os.listdir(app_root)

                if 'precomputed' in app_root_contents:
                    output.debug.set(f"precomputed/ FOUND in {app_root}")
                    precomputed_path = os.path.join(app_root, 'precomputed')
                    output.debug.set(f"Files in precomputed/: {os.listdir(precomputed_path)}")
                else:
                    output.debug.set(f"precomputed/ NOT found in {app_root}. Contents: {app_root_contents}")
                diff_folder, new_diff = find_closest_diff_folder(n_folder_path, differentiate)
                output_text_1=f'Something went wrong while looking for your file \n'
                output_text_2=f'Use the command below to calculate a pooling strategy for {n_samp} samples with up to {differentiate} positives'
                output_text=output_text_1+output_text_2
                output.database_reply.set(output_text)
                output.summary_table.set(pd.DataFrame(columns=ls_met))
                output.allow_downloads.set(False)
                
            
            
        # Prepare the correct command for pool_N.py based on its arguments
        command_p = f"python pool_N.py --n_samp {n_samp} --differentiate {differentiate} --method all --path your/path"
        intro_command='You can always run the code locally for your specific case, for the values you selected you need to run:\n\n'
        output.personalized_command.set(intro_command+command_p)

        decode_p = f"python decode_N.py --differentiate {differentiate} path_to_WA your/path/to/pooling_strategy.csv --readout path/to/readout.csv (or your readout e.g. 1,0,0,1,0)"
        intro_command='You can must run the decoder locally for your specific pooling strategy and readout:\n\n'
        output.decoder.set(intro_command+decode_p)

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
                def debug():
                    return  output.debug.get()
                



    @reactive.Effect
    @reactive.event(input.prev_submit)
    def _():
        n = input.prev_n_samp()
        p = input.prev_prevalence()
        max_e = input.prev_max_error()
        output.prev_n_samp.set(n)
        output.prev_prevalence.set(p)
        output.prev_max_error.set(max_e)
        output.last_val_prev.set(f"Number of Samples: {n}, Prevalence: {p}, Max error: {max_e}")

    @output
    @render.text
    def last_val_prev():
        return output.last_val_prev.get()


    @output
    @render.text
    def allow_downloads():
        return str(output.allow_downloads.get()).lower()  # returns 'true' or 'false'

    @output
    @render.text
    def last_val():
        return  output.last_values.get()
    
    @output
    @render.text
    def database_r():
        return  output.database_reply.get()
    
    @output
    @render.ui
    def summary_text():
        return ui.HTML(output.summary_text.get())
    
    @output
    @render.text
    def commands():
        return  output.personalized_command.get()
    
    @output
    @render.text
    def decode():
        return  output.decoder.get()
    
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
        yield DFFS['Matrix'].to_csv(index=True)
    
    @output
    @render.download(filename=lambda: "2D_pooling.csv")
    async def download_table_2d():
        # Yield the content of the CSV file
        DFFS=output.dataframes.get()
        #idt=[i for i in list(DFFS) if i.startswith('multidim')]
        yield DFFS['multidim-2'].to_csv(index=True)

    @output
    @render.download(filename=lambda: "3D_pooling.csv")
    async def download_table_3d():
        # Yield the content of the CSV file
        DFFS=output.dataframes.get()
        #idt=[i for i in list(DFFS) if i.startswith('multidim')]
        yield DFFS['multidim-3'].to_csv(index=True)

    @output
    @render.download(filename=lambda: "4D_pooling.csv")
    async def download_table_4d():
        # Yield the content of the CSV file
        DFFS=output.dataframes.get()
        #idt=[i for i in list(DFFS) if i.startswith('multidim')]
        yield DFFS['multidim-4'].to_csv(index=True)

    @output
    @render.download(filename=lambda: "random_pooling.csv")
    async def download_table_random():
        # Yield the content of the CSV file
        DFFS=output.dataframes.get()
        yield DFFS['Random'].to_csv(index=True)
    
    @output
    @render.download(filename=lambda: "STD_pooling.csv")
    async def download_table_STD():
        # Yield the content of the CSV file
        DFFS=output.dataframes.get()
        yield DFFS['STD'].to_csv(index=True)

    @output
    @render.download(filename=lambda: "Chinese_remainder_pooling.csv")
    async def download_table_CT():
        # Yield the content of the CSV file
        DFFS=output.dataframes.get()
        yield DFFS['Chinese trick'].to_csv(index=True)

    @output
    @render.download(filename=lambda: "Chinese_bktrk_pooling.csv")
    async def download_table_CTB():
        # Yield the content of the CSV file
        DFFS=output.dataframes.get()
        yield DFFS['Chinese bktrk'].to_csv(index=True)

    @output
    @render.download(filename=lambda: "Chinese_special_pooling.csv")
    async def download_table_CTS():
        # Yield the content of the CSV file
        DFFS=output.dataframes.get()
        if 'Chinese special' in DFFS.keys():
            yield DFFS['Chinese special'].to_csv(index=True)
        else:
            yield DFFS['Chinese bktrk'].to_csv(index=True)

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
    @render.ui
    def prevalence_error_table():
        N = output.prev_n_samp.get()
        prevalence = output.prev_prevalence.get()
        max_err = output.prev_max_error.get()
        df_orig = expected_error_table(N, prevalence, max_diff=4, correct=False)
        df_corr = expected_error_table(N, prevalence, max_diff=4, correct=True)

        def color_for_value(val_orig, val_corr):
            if val_orig > max_err and val_corr > max_err:
                return "background-color: #ffeaea; color: #111;"  # red
            elif val_corr > max_err:
                return "background-color: #ffffea; color: #111;"  # yellow
            else:
                return "background-color: #b6fcb6; color: #111;"  # green

        # Build two tables side by side
        html = '<div style="display: flex; flex-direction: row; gap: 32px; justify-content: center;">'

        # Original table
        html += '<div>'
        html += '<div style="text-align: center; font-weight: bold; margin-bottom: 8px;">Single batch</div>'
        html += '<table style="border-collapse: collapse; font-size: 18px; min-width: 300px;">'
        html += '<tr><th style="border: 1px solid #aaa; padding: 10px 16px; text-align: center;">N \\ Differentiate</th>'
        for col in df_orig.columns:
            html += f'<th style="border: 1px solid #aaa; padding: 10px 16px; text-align: center;">{col}</th>'
        html += '</tr>'
        for idx, row in df_orig.iterrows():
            html += f'<tr><th style="border: 1px solid #aaa; padding: 10px 16px; text-align: center;">{idx}</th>'
            for j, val in enumerate(row):
                val_corr = df_corr.loc[idx, df_corr.columns[j]] if idx in df_corr.index and df_corr.columns[j] in df_corr.columns else None
                style = color_for_value(val, val_corr)
                html += f'<td style="border: 1px solid #aaa; padding: 10px 16px; text-align: center; {style}">{val:.3g}</td>'
            html += '</tr>'
        html += '</table>'
        html += '</div>'

        # Corrected table
        html += '<div>'
        html += '<div style="text-align: center; font-weight: bold; margin-bottom: 8px;">Corrected for FWER</div>'
        html += '<table style="border-collapse: collapse; font-size: 18px; min-width: 300px;">'
        html += '<tr><th style="border: 1px solid #aaa; padding: 10px 16px; text-align: center;">N \\ Differentiate</th>'
        for col in df_corr.columns:
            html += f'<th style="border: 1px solid #aaa; padding: 10px 16px; text-align: center;">{col}</th>'
        html += '</tr>'
        for idx, row in df_corr.iterrows():
            html += f'<tr><th style="border: 1px solid #aaa; padding: 10px 16px; text-align: center;">{idx}</th>'
            for j, val_corr in enumerate(row):
                val_orig = df_orig.loc[idx, df_orig.columns[j]] if idx in df_orig.index and df_orig.columns[j] in df_orig.columns else None
                style = color_for_value(val_orig, val_corr)
                html += f'<td style="border: 1px solid #aaa; padding: 10px 16px; text-align: center; {style}">{val_corr:.3g}</td>'
            html += '</tr>'
        html += '</table>'
        html += '</div>'

        html += '</div>'  # end flex
        return ui.HTML(html)

    @output
    @render.ui
    def prevalence_legend():
        html = '<div style="font-size: 15px; text-align: center; margin-top: 18px;">'
        html += '<div style="margin-bottom: 6px;"><span style="background-color: #ffeaea; color: #111; padding: 4px 14px; border-radius: 4px; border: 1px solid #eee;">Red: pooling error &gt; max error for both single pooling and FWER (all poolings combined)</span></div>'
        html += '<div style="margin-bottom: 6px;"><span style="background-color: #ffffea; color: #111; padding: 4px 14px; border-radius: 4px; border: 1px solid #eee;">Yellow: pooling error ≤ max error for single pooling, but &gt; max error for FWER (all poolings combined)</span></div>'
        html += '<div><span style="background-color: #b6fcb6; color: #111; padding: 4px 14px; border-radius: 4px; border: 1px solid #eee;">Green: pooling error ≤ max error for both single pooling and FWER (all poolings combined)</span></div>'
        html += '</div>'
        return ui.HTML(html)

    @output
    @render.ui
    def prevalence_explanation():
        html = '<div style="margin-top: 28px; font-size: 14px; color: #444; text-align: center;">These tables are a guide to help you decide what kind of pooling strategies might be right for your particular problem.</div>'
        return ui.HTML(html)



'''
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
'''

# Create the app object
app = App(app_ui, server)