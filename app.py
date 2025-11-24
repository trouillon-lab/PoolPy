from shiny import App, ui, reactive, render
from pathlib import Path
import pandas as pd
import os
import pickle
import numpy as np
import os
import pandas as pd
import scipy.stats
import re
import itertools
import string
#import matplotlib.pyplot as plt
summary_text = ""
# Define the UI
app_ui = ui.page_fluid(
    ui.tags.title("PoolPy"),
    ui.navset_tab(
        ui.nav_panel("Pooling Methods Comparison",
            # App logo/image at the top
            ui.div(
                ui.output_image("logo_img_main", height='0px'),
                style="display: block; padding: 0; margin: 0; text-align: centre;"
            ),
            # Centered Title
            ui.h2("Pooling Methods Comparison", style="text-align: center; margin-top: 150px;"),
            # Short description below the title
            ui.div(
                "This app evaluates the performances of different pooling methods. Designs for specific use cases can be downloaded below.",
                ui.br(),
                "For details, see the ",
                ui.a("Methods and Metrics", href="#methods-metrics", id="nav-methods-link", style="text-decoration: underline; cursor: pointer;"),
                " section.",
                " For information on the code, visit our ",
                ui.a("GitHub", href="https://github.com/trouillon-lab/PoolPy", style="text-decoration: underline;", target="_blank", rel="noopener noreferrer"),
                ".",
                ui.br(),
                "To determine a max. number of positives based on prevalence, consult the ",
                ui.a("Prevalence", href="#prevalence-section", id="nav-prevalence-link", style="text-decoration: underline; cursor: pointer;"),
                " section.",
                ui.br(),
                "To decode results from a pooled experiment, consult the ",
                ui.a("Decoder", href="#decoder-section", id="nav-decoder-link", style="text-decoration: underline; cursor: pointer;"),
                " section.",
                ui.br(),
                "To generate a method file for a pipetting robot following a pooling design, consult the ",
                ui.a("Automation", href="#automation-section", id="nav-automation-link", style="text-decoration: underline; cursor: pointer;"),
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
                    var decoderLink = document.getElementById('nav-decoder-link');
                    if (decoderLink) {
                        decoderLink.addEventListener('click', function(event) {
                            event.preventDefault();
                            var tab = document.querySelector("[data-value='Decoder']");
                            if (tab) tab.click();
                        });
                    }
                    var automationLink = document.getElementById('nav-automation-link');
                    if (automationLink) {
                        automationLink.addEventListener('click', function(event) {
                            event.preventDefault();
                            var tab = document.querySelector("[data-value='Automation']");
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
                ui.column(2, ui.input_numeric("differentiate", "Max. number of positives:", value=1)),
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
                ui.output_ui("last_val"),
                style="text-align: center;"
            ),
            ui.hr(),
            ui.div(
                ui.h4("Database reply:"),
                ui.output_ui("database_r"),
                style="text-align: center;"
            ),
            ui.hr(),
            ui.panel_conditional(
                "output.allow_downloads == 'true'",
                ui.div(
                    ui.div(
                        ui.output_ui("summary_table_title"),
                        style="margin: 32px 0 32px 0;"
                    ),
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
                        style="text-align: center; margin-top: 20px;"
                    )
                )
            ),
            ui.hr(),
            ui.div(
                ui.h4("Downloads & code"),
                style="text-align: center;"
            ),
            ui.div("", style="height: 32px;"),
            ui.panel_conditional(
                "output.allow_fly == 'true'",
                ui.div(
                    ui.output_ui("fly_download_dip_title"),
                    ui.download_button("download_table_STD_fly", "STD"),
                    ui.download_button("download_table_CT_fly", "Chinese Remainder"),
                    ui.download_button("download_table_CT_bktrk_fly", "Ch. Rm. Backtrack"),
                    ui.download_button("download_table_CT_special_fly", "Ch. Rm. Special"),
                    style="text-align: center;",
                ),
                ui.div("", style="height: 20px;")
            ),
            ui.div("", style="height: 32px;"),
            ui.panel_conditional(
                "output.allow_fly == 'true'",
                ui.div(
                    ui.output_ui("fly_download_ind_title"),
                    ui.download_button("download_table_matrix_fly", "Matrix"),
                    ui.download_button("download_table_2d_fly", "2-dim"),
                    ui.download_button("download_table_3d_fly", "3-dim"),
                    ui.download_button("download_table_4d_fly", "4-dim"),
                    ui.download_button("download_table_binary_fly", "Binary"),
                    style="text-align: center;",
                ),
                ui.div("", style="height: 20px;")
            ),
            ui.panel_conditional(
                "output.allow_downloads == 'true'",
                ui.div(
                    ui.h4("Downloadable precomputed random design", style="margin-bottom: 20px;"),
                    #ui.download_button("download_table_matrix", "Matrix"),
                    #ui.download_button("download_table_2d", "2-dim"),
                    #ui.download_button("download_table_3d", "3-dim"),
                    #ui.download_button("download_table_4d", "4-dim"),
                    ui.download_button("download_table_random", "Random"),
                    #ui.download_button("download_table_STD", "STD"),
                    #ui.download_button("download_table_CT", "Chinese remainder"),
                    #ui.download_button("download_table_CTB", "Chinese rm bktrk"),
                    #ui.download_button("download_table_CTS", "Chinese rm special"),
                    #ui.download_button("download_table_binary", "Binary"),
                    style="text-align: center;",
                ),
                ui.div("", style="height: 20px;"),
            ),
            ui.div(
                ui.h4("Command to generate designs locally"),
                ui.output_text_verbatim("commands"),
                style="text-align: center;"
            ),
            ui.panel_conditional(
                "output.allow_fly == 'true'",
                ui.div(
                    ui.h4("Command to decode results locally"),
                    ui.output_text_verbatim("decode"),
                    style="text-align: center;"
                ),
            ),
            ui.div(
                ui.output_text_verbatim("allow_downloads"),
                style="position: fixed; bottom: 2px; right: 2px; opacity: 0.05; font-size: 8px; z-index: 1; pointer-events: none;"
            ),
            ui.div(
                ui.output_text_verbatim("allow_fly"),
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
        ui.nav_panel("Decoder",
            ui.div(
                ui.output_image("logo_img_decoder", height='0px'),
                style="display: block; padding: 0; margin: 0; text-align: centre;"
            ),
            ui.h2("Decoder", style="text-align: center; margin-top: 150px;"),
            
            # Short description below the title
            ui.div(
                "Results from tests done following a pooling design from PoolPy can be decoded here.",
                ui.br(),
                "First, load your design file. Only csv designs files generated by PoolPy are accepted.",
                ui.br(),
                "Then, enter in the readout which pools were positive in your tests as a list of positive pools (e.g. 1,3,5).",
                style="text-align: center; margin-bottom: 18px; font-size: 16px; color: #444;"
            ),
            
            ui.div(
                ui.input_file("uploaded_csv_decoder", "Upload your design csv file (max. 100kb)", accept=[".csv"], multiple=False),
                style="display: flex; justify-content: center; margin-bottom: 22px;"
            ),
            ui.div(
                ui.input_text(
                    "readout_string",
                    "Readout (positive pools):",
                    value="",
                    placeholder="Enter comma-separated values"
                ),
                style="display: flex; justify-content: center; margin-bottom: 18px;"
            ),
            ui.div(
                ui.input_text(
                    "decoder_diff",
                    "Max. number of positives (optional):",
                    placeholder="Leave empty for auto inference"
                ),
                style="display: flex; justify-content: center; margin-bottom: 18px;"
            ),
            ui.div(
                ui.input_action_button("readout_submit", "Enter", style="display: block; margin: 0 auto; margin-bottom: 32px;"),
                style="text-align: center;"
            ),
            ui.tags.script('''
                document.addEventListener('DOMContentLoaded', function() {
                    ["readout_string", "decoder_diff"].forEach(function(id) {
                        var el = document.getElementById(id);
                        if (el) {
                            el.addEventListener('keydown', function(event) {
                                if (event.key === 'Enter') {
                                    event.preventDefault();
                                    var btn = document.getElementById('readout_submit');
                                    if (btn) btn.click();
                                }
                            });
                        }
                    });
                });
            '''),
            ui.div(
                ui.output_ui("database_r_decoder"),
                style="text-align: center; margin-bottom: 18px; font-size: 16px; color: #444;"
            ),
            ui.panel_conditional(
                "output.allow_decoder == 'true'",
                ui.div(
                    ui.download_button("download_decoder_output", "Download decoder output"),
                    style="display: flex; justify-content: center; gap: 32px; margin-bottom: 32px;"
                ),
            ),
            ui.div(
                ui.output_text_verbatim("allow_decoder"),
                style="position: fixed; bottom: 2px; right: 2px; opacity: 0.05; font-size: 8px; z-index: 1; pointer-events: none;"
            ),
        ),

        ui.nav_panel("Methods & Metrics",
            ui.div(
                ui.output_image("logo_img_methods", height='0px'),
                style="display: block; padding: 0; margin: 0; text-align: cenrte;"
            ),
            ui.h2("Methods & Metrics",  style="text-align: center; margin-top: 150px;"),
            ui.div(
                ui.h4("Differentiate-dependent Pooling Methods:"),
                ui.div(
                    "These methods change their design also based on the number of expected positive samples. They often are non-adaptive strategy and are well suited for high differentiate values.",
                    style="margin-bottom: 10px; font-size: 15px; color: #444;"
                ),
                ui.tags.ul(
                    ui.tags.li(ui.tags.b("Random"), ": Semi-adaptive method where pools are formed by randomly assigning samples."),
                    ui.tags.li(ui.tags.b("STD"), ": Non-adaptive method based on prime numbers and modulus operations."),
                    ui.tags.li(ui.tags.b("Chinese Remainder methods"), ": Non-adaptive method based on the Chinese Remainder Theorem, with variants for backtracking (backtrack) and special cases (special)."),
                    ui.tags.li(ui.tags.b("Hierachical"), ": Strictly adaptive method that iteratively partitions the set of possible positive samples." \
                    "For this method, the N pools metric lists the number of splits at each stage. " \
                    "For example, [3,3] means first divide samples into 3 pools (as equal as possible), then test and split any positive pool into 3 again." \
                    "After testing of these pools, each sample from each positive pool finally needs to be individually tested." \
                    ),
                ),

                ui.h4("Differentiate-independent Pooling Methods:"),
                ui.div(
                    "These methods do not change their output design based on the number of expected positive samples. They are to be used with caution if more than 1 positive sample is expected.",
                    style="margin-bottom: 10px; font-size: 15px; color: #444;"
                ),
                ui.tags.ul(
                    ui.tags.li(ui.tags.b("Matrix"), ": Semi-adaptive method where each sample is included in a unique row and column pool."),
                    ui.tags.li(ui.tags.b("Multidimensional (3D, 4D)"), ": Semi-adaptive method where samples are arranged in higher-dimensional matrices, each coordinate in each dimension representing a pool."),
                    ui.tags.li(ui.tags.b("Binary"), ": Semi-adaptive method where samples are assigned to pools according to a binary code, maximizing information per test."),
                ),
                ui.h4("Metrics in the Summary Table:"),
                ui.tags.ul(
                    ui.tags.li(ui.tags.b("Method"), ": Name of the pooling method used."),
                    ui.tags.li(ui.tags.b("Mean experiments"), ": Average number of tests required to identify the positive samples."),
                    ui.tags.li(ui.tags.b("Max. samples per pool"), ": Maximum number of samples combined in any single pool."),
                    ui.tags.li(ui.tags.b("N pools"), ": Total number of pools used in the first step of the strategy. For Hierachical this is a list of splits, explained above."),
                    ui.tags.li(ui.tags.b("Percentage check"), ": Fraction of cases requiring additional verification or retesting beyond the first step."),
                    ui.tags.li(ui.tags.b("Mean extra experiments"), ": Average number of extra tests needed beyond the first step.")
                ),
                ui.h4("Notation:"),
                ui.tags.ul(
                    ui.tags.li(ui.tags.b("D"), ": Differentiate, maximum number of positive samples."),
                    ui.tags.li(ui.tags.b("S"), ": Number of samples to test."),
                    ui.tags.li(ui.tags.b("ρ"), ": Prevalence of positives in the population."),
                    ui.tags.li(ui.tags.b("W"), ": Total number of pools needed for a method."),
                ),
                style="max-width: 700px; margin: 0 auto; font-size: 16px; color: #333;"
            )
        ),
        ui.nav_panel("Prevalence",
                ui.div(
                    ui.output_image("logo_img_prev", height='0px'),
                    style="display: block; padding: 0; margin: 0; text-align: centre;"
                ),
                ui.h2("Prevalence", style="text-align: center; margin-top: 150px;"),
                ui.div(
                    ui.div(
                        "This section is meant to guide the decision of pooling parameters for specific use cases. Particularly, this is to choose a max. number of positive samples and a number of test batches.",
                        ui.br(),
                        "The two tables reflect cases where all samples are processsed into one single batch (left) or broken down into multiple batches (right).",
                        ui.br(),
                        "The tables show the probability of making at least one mistake while reading the results of either one (left) or multiple combined (right) combinatorial pooling batch(es).",
                        ui.br(),
                        "Based on a provided prevalence estimate, error rates are reported over ranges of sample (and batch) numbers (S; rows) and of max. number of positive samples (D; columns).",
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
                    ui.output_ui("last_val_prev"),
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
                #redundant?
                # ui.div(
                #     ui.output_ui("prevalence_explanation"),
                #     style="margin-top: 8px; display: flex; justify-content: center;"
                # ),
            ),
        
        ui.nav_panel("Automation",
            ui.div(
                ui.output_image("logo_img_automation", height='0px'),
                style="display: block; padding: 0; margin: 0; text-align: centre;"
            ),
            ui.h2("Automation", style="text-align: center; margin-top: 150px;"),
            
            # Short description below the title
            ui.div(
                "This tool generates a method file for a pipetting robot to follow a pooling design file.",
                ui.br(),
                "Load your design file below. Only csv designs files generated by PoolPy are accepted.",
                ui.br(),
                "The generated file lists all pipetting steps required to generate the pools following the vendor's format. If your brand is not supported, reach out to us.",
                ui.br(),
                "Volumes and plate names might have to be adjusted to your specific use case.",
                style="text-align: center; margin-bottom: 18px; font-size: 16px; color: #444;"
            ),
            
            ui.div(
                ui.input_file("uploaded_csv_auto", "Upload your pooling strategy csv file (max. 100kb)", accept=[".csv"], multiple=False),
                style="display: flex; justify-content: center; margin-bottom: 32px;"
            ),
            ui.div(
                ui.output_ui("database_r_auto"),
                style="text-align: center; margin-bottom: 18px; font-size: 16px; color: #444;"
            ),
            ui.panel_conditional(
                "output.allow_robot == 'true'",
                ui.div(
                    ui.download_button("download_biomek", "Download Biomek robot code"),
                    ui.download_button("download_hamilton", "Download Hamilton robot code"),
                    ui.download_button("download_opentron", "Download Opentron robot code"),
                    style="display: flex; justify-content: center; gap: 32px; margin-bottom: 32px;"
                ),
            ),
            ui.div(
                ui.output_text_verbatim("allow_robot"),
                style="position: fixed; bottom: 2px; right: 2px; opacity: 0.05; font-size: 8px; z-index: 1; pointer-events: none;"
            ),
        ),

    )
)



WA_SUB_DIRECTORY='precomputed'
#SCRAMBLER_DIRECTORY='.\output'
MAX_DIFFERENTIATE=4
MAX_N=1000


def from_id_to_well(id, offset=0):
    plate=id//96
    id=id-plate*96
    row=id//12
    column=id-row*12
    well=chr(65 + row)
    return str(plate+1+offset), str(well)+str(column+1)


def clean_WA(b):
    b1=b.astype(int)
    tmp1=pd.DataFrame(b1, columns=['Pool '+ str(i) for i in range(b1.shape[1])], index=['Sample '+ str(i) for i in range(b1.shape[0])])
    return tmp1


def int_to_base(n, N):
    """ Return base N representation for int n. """
    base_n_digits = string.digits + string.ascii_lowercase + string.ascii_uppercase
    result = ""
    if n < 0:
        sign = "-"
        n = -n
    else:
        sign = ""
    while n > 0:
        q, r = divmod(n, N)
        result += base_n_digits[r]
        n = q
    if result == "":
        result = "0"
    return sign + "".join(reversed(result))


#Coumpound counter starts from 1
# Helper function for binary translation
def IntegerToBinaryTF(num: int, ls_bn: list)-> list:
    if num >= 2:
        ls_bn=IntegerToBinaryTF(num // 2, ls_bn)
    ls_bn.append(num % 2==1)
    return(ls_bn)

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))

def isprime(n:int):
    return re.compile(r'^1?$|^(11+)\1+$').match('1' * n) is None

def get_Gamma(q,N):
    return(int(np.ceil(np.log(N)/np.log(q))-1))

def get_s(N,j,q):
    vec=np.arange(N)
    out_vec=vec.copy()
    Gamma=get_Gamma(q,N)
    for ct in range(Gamma):
        c=ct+1
        out_vec=out_vec+j**c*(vec//q**c)
    return(out_vec)

def STD(N,q,k):
    L=np.zeros((k,q,N))==1
    for j in range(k):
        s=get_s(N,j,q)
        s=s%q
        idc=np.arange(N)
        L[j,s,idc]=True

    L=L.reshape(-1,N)
    return(L.T)


# Function to assign each compound to its wells
def well_selecter(compound: int, n_wells:int, differentiate=1) -> np.array:
    if differentiate not in [1,2]:
        print('For the moment this code is only able to create well assignments matrices to distinguish up to combinations of 2 active compounds')
        return(-1)
    if differentiate==1:
        ls_bn=[]
        used_wells=IntegerToBinaryTF(compound, ls_bn)
        sel_wells=[False]*(n_wells-len(used_wells))+used_wells
    if differentiate==2:
        for i in range((n_wells-1)//3+1):
            if 0<compound and compound <= n_wells-1-3*i:
                sel_wells=(n_wells-1-3*i-compound)*[False]+[True]+i*3*[False]+[True]+[False]*(compound-1)
                break
            compound=compound-(n_wells-1-3*i)
    return(np.array(sel_wells))

def assign_wells_mat(n_compounds:int, **kwargs)->np.array:
    L1=np.ceil(np.sqrt(n_compounds))
    L2=L1-1 if L1*(L1-1)>=n_compounds else L1
    well_assigner=np.zeros((n_compounds, int(L1+L2)))==1
    for i in range(n_compounds):
        cp_id=[int(i//L2), int(L1+(i % L2))]
        well_assigner[i,cp_id]=True
    return(well_assigner)

# This functions also identifies the minimum number of wells needed for the compounds and level of detail (differentiate) selected
def assign_wells_bin(n_compounds:int, differentiate=1, **kwargs) -> np.array:

    n_wells=int(np.ceil(np.log2(n_compounds +1)))
        
    well_assigner=np.zeros((n_compounds, n_wells))==1
    for i in range(n_compounds):
        well_assigner[i,:]=well_selecter(i+1, n_wells, differentiate=1)
    return(well_assigner)


''' Method 3: Pooling using multidimensional matrix design'''
def assign_wells_multidim(n_compounds:int, n_dims:int, **kwargs)->np.array:
    L1=np.ceil(np.power(n_compounds, 1/n_dims))
    i=0
    while np.power(L1, n_dims-i-1)*np.power(L1-1, i+1)>=n_compounds:
        i=i+1
    ls_dim=[L1]*(n_dims-i)+[L1-1]*i
    up_samps=np.prod(np.array(ls_dim))
    well_assigner=np.zeros((n_compounds, int(L1*(n_dims-i)+(L1-1)*i)))==1
    for j in range(n_compounds):
        cp_id=[]
        jj=np.copy(j)
        rem_dim=up_samps
        past_dims=0
        for k in range(n_dims):
            rem_dim=rem_dim/ls_dim[k]
            js=jj//rem_dim
            jj=jj-js*rem_dim
            jd=js+past_dims
            cp_id.append(int(jd))
            past_dims=past_dims+ls_dim[k]
        well_assigner[j,cp_id]=True
    return(well_assigner)

def assign_wells_STD(n_compounds:int, differentiate=1, False_results=0, force_q=False, **kwargs):
    N=n_compounds
    t=differentiate
    E=False_results
    poss_q=[x for x in range(n_compounds) if isprime(x)]
    for q in poss_q:
        if t*get_Gamma(q,N)+2*E<q:
            break
    if isprime(force_q):
        if t*get_Gamma(force_q,N)+2*E<force_q:
            q=force_q 
    Gamma=get_Gamma(q,N)
    k=t*Gamma+2*E+1
    WA=STD(N,q,k)
    return(WA)

def assign_wells_chinese_old(n_compounds:int,  differentiate:int,**kwargs)->np.array:
    prod=1
    n=1
    primes=[]
    c_id=np.arange(n_compounds) 
    while prod<n_compounds**differentiate:
        n=n+1
        if isprime(n):
            prod=prod*n
            primes.append(n)
    WA=np.zeros((np.sum(primes), n_compounds))==1
    past_primes=0
    for prime in primes:
        temp_wa=np.zeros((prime, n_compounds))==1
        for x in range(prime):
            ids=c_id%prime==x    
            temp_wa[x, ids]=True
        WA[past_primes:past_primes+prime,:]=temp_wa
        past_primes=past_primes+prime

    return(WA.T)


def assign_wells_chinese(n_compounds:int,  differentiate:int, backtrack=False, special_diff=False, **kwargs)->np.array:
    prod=1
    n=1
    primes=[]
    c_id=np.arange(n_compounds) 
    while prod<n_compounds**differentiate:
        n=n+1
        if isprime(n):
            prod=prod*n
            primes.append(n)

    if backtrack:
        T=np.inf
        nprimes=np.array(primes)
        ND=n_compounds**differentiate
        ls_of_ls=[]
        LMP=np.log(primes[-1])
        for pi in primes:
            LE=np.floor(LMP/np.log(pi)).astype(int)
            ls_of_ls.append(list(range(LE+1)))
        ls_iter=list(itertools.product(*ls_of_ls))
        for id_combo, combo in enumerate(ls_iter):
            carr=np.array(combo)
            flt=carr>0
            this_primes=nprimes[flt]
            this_exp=carr[flt]
            npc=this_primes**this_exp
            if np.prod(npc)>=ND and np.sum(npc)<T:
                T=np.sum(npc)
                best_id=id_combo
        combo=ls_iter[best_id]
        carr=np.array(combo)
        flt=carr>0
        this_primes=nprimes[flt]
        this_exp=carr[flt]
        npc=this_primes**this_exp

        WA=np.zeros((np.sum(npc), n_compounds))==1
        past_primes=0
        for prime in npc:
            temp_wa=np.zeros((prime, n_compounds))==1
            for x in range(prime):
                ids=c_id%prime==x    
                temp_wa[x, ids]=True
            WA[past_primes:past_primes+prime,:]=temp_wa
            past_primes=past_primes+prime

        return(WA.T)   
    
    if special_diff and differentiate==2:
        q=np.ceil(np.log(n_compounds)/np.log(3)).astype(int)
        t=int((q+5)*q/2)
        WA=np.zeros((t, n_compounds))==1
        ls_nc3=[list(i)[::-1] for i in [int_to_base(j,3).zfill(q) for j in range(n_compounds)]]
        for i in range(q):
            for ii in range(3):
                for j in range(n_compounds):
                    WA[3*i+ii,j]=True if int(ls_nc3[j][i])==int(ii) else False
        k=3*q
        for i in range(q):
            for ii in range(i+1,q):
                for j in range(n_compounds):
                    WA[k,j]=True if int(ls_nc3[j][i])==int(ls_nc3[j][ii]) else False
                k+=1
        return(WA.T)
    
    if special_diff and differentiate==3:
        q=np.ceil(np.log(n_compounds)/np.log(2)).astype(int)
        t=int((q-1)*q*2)
        WA=np.zeros((t, n_compounds))==1
        ls_nc3=[list(i)[::-1] for i in [int_to_base(j,2).zfill(q) for j in range(n_compounds)]]
        k=0
        for i in range(q):
            for ii in range(i+1,q):
                for nu in [0,1]:
                    for nuu in [0,1]:
                        for j in range(n_compounds):
                            WA[k,j]=True if int(ls_nc3[j][i])==nu and int(ls_nc3[j][ii])==nuu else False
                        k+=1
        return(WA.T)


    WA=np.zeros((np.sum(primes), n_compounds))==1
    past_primes=0
    for prime in primes:
        temp_wa=np.zeros((prime, n_compounds))==1
        for x in range(prime):
            ids=c_id%prime==x    
            temp_wa[x, ids]=True
        WA[past_primes:past_primes+prime,:]=temp_wa
        past_primes=past_primes+prime

    return(WA.T)
 
def decode_precomp(well_assigner:np.array, differentiate:int, 
                   scrambler:dict, readout:np.ndarray, **kwargs) -> list:
    if differentiate==0:
        return(True,well_assigner, np.array([1]*well_assigner.shape[0]))
    
    N=well_assigner.shape[0]
    sc_list=np.arange(N).tolist()
    for i in range(differentiate):
        diff=i+1
        if diff ==1:
            full_well_assigner=well_assigner.copy()
        else:
            this_sc=scrambler[diff]
            full_well_assigner=np.concatenate((full_well_assigner,np.any(well_assigner[this_sc], axis=1)))
            sc_list.extend(this_sc.tolist())


    idxs = np.all(readout == full_well_assigner, axis=1)
    return list(itertools.compress(sc_list,idxs))
        
    

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
    DF=None
    while DF is None:
        DF=find_closest_diff_folder_inner(n_folder_path=n_folder_path, differentiate=differentiate)
        differentiate=differentiate-1
    return DF




def find_closest_diff_folder_inner(n_folder_path, differentiate):
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


def find_closest_n_diff_folder(wa_directory, n_samp, diff):
    # Step 1: Find all N_X folders
    n_folders = [f for f in os.listdir(wa_directory) if os.path.isdir(os.path.join(wa_directory, f)) and f.startswith('N_')]
    n_x_values = [int(f.split('_')[1]) for f in n_folders]
    n_folder_map = {int(f.split('_')[1]): f for f in n_folders}

    # Step 2: For each N_X, find all diff_Y subfolders and collect all available differentiate values
    #diff_values = [[] for i in range(MAX_DIFFERENTIATE)]
    diff_values = set()
    diff_map = {}  # (diff, N) -> folder path
    for n_val, n_folder in n_folder_map.items():
        n_folder_path = os.path.join(wa_directory, n_folder)
        diff_folders = [f for f in os.listdir(n_folder_path) if os.path.isdir(os.path.join(n_folder_path, f)) and f.startswith('diff_')]
        for d_folder in diff_folders:
            try:
                y = int(d_folder.split('_')[1])
                diff_folder_path = os.path.join(n_folder_path, d_folder)
                # Check for a file starting with 'Metrics' in this diff folder
                metrics_files = [mf for mf in os.listdir(diff_folder_path) if mf.startswith('Metrics')]
                if metrics_files:
                    diff_values.add(y)
                    diff_map[(y, n_val)] = diff_folder_path
            except Exception:
                continue

    # Step 3: Find the best differentiate value to use
    available_diffs = sorted(diff_values)
    if not available_diffs:
        return None  # No differentiate folders found at all

    if diff in available_diffs:
        chosen_diff = diff
    else:
        less_than = [d for d in available_diffs if d < diff]
        if less_than:
            chosen_diff = max(less_than)
        else:
            chosen_diff = max(available_diffs)

    # Step 4: For chosen_diff, find all N_X with that diff
    n_with_diff = [n for (d, n) in diff_map if d == chosen_diff]
    if not n_with_diff:
        return None  # No folder found

    # Step 5: Find the closest N >= n_samp, else closest overall
    n_with_diff_sorted = sorted(n_with_diff)
    n_candidates = [n for n in n_with_diff_sorted if n >= n_samp]
    if n_candidates:
        chosen_n = min(n_candidates)
    else:
        # If none >= n_samp, take the closest overall
        chosen_n = min(n_with_diff_sorted, key=lambda x: abs(x - n_samp))

    # Step 6: Return the folder path and the chosen values
    folder_path = diff_map[(chosen_diff, chosen_n)]
    return folder_path, chosen_n, chosen_diff


# Function to generate expected error table for prevalence tab
def expected_error_table(N, prevalence, max_diff=None, correct=False, max_error=0.05, extra_steps=2):
    if max_diff is None:
        max_diff = 4  # or set to N or another reasonable value
    diff_values = np.arange(1, int(max_diff)+extra_steps+1)
    extreme_diff= int(np.ceil(scipy.stats.binom.isf(max_error, N, prevalence)))
    min_N=5
    if extreme_diff>max_diff+extra_steps:
        rooty=np.power(extreme_diff/max_diff,1/(extra_steps))
        for i in range(extra_steps):
            new_diff=int(np.ceil(max_diff*np.power(rooty, i+1)))
            diff_values[max_diff+i]=new_diff
    #elif max_diff<extreme_diff<max_diff+extra_steps:
    #    diff_values = np.arange(1, int(extreme_diff)+1)
    else:
        diff_values = np.arange(1, int(np.max([max_diff, extreme_diff]))+1)

    pop_sizes = []
    N_current = int(np.ceil(N))
    factor=[]
    N_factor=[]
    NF=1
    # Build the population sizes as described, always take the ceiling to be on the safe side
    while N_current >= min_N:
        pop_sizes.append(N_current)
        factor.append(NF)
        N_factor.append(str(NF)+' X '+str(N_current))
        N_current = int(np.ceil(N_current / 2))
        NF=int(NF*2)
        if N_current < min_N:
            break
        N_factor.append(str(NF)+' X '+str(N_current))
        pop_sizes.append(N_current)
        factor.append(NF)
        NF=int(NF*2.5)
        N_current = int(np.ceil(N_current / 2.5))
        if N_current < min_N:
            break
        N_factor.append(str(NF)+' X '+str(N_current))
        pop_sizes.append(N_current)
        factor.append(NF)
        NF=int(NF*2)
        N_current = int(np.ceil(N_current / 2))
    # Remove duplicates and sort descending
    #pop_sizes = sorted(set([x for x in pop_sizes if x >= min_N]), reverse=True)
    
    # Build the table
    
    data = []
    CF=pop_sizes
    for N_val,NF in zip(pop_sizes,factor):
        row = []
        for diff in diff_values:
            p_error = scipy.stats.binom.sf(diff, N_val, prevalence)
            if correct:
                p_error=1-(1-p_error)**int(NF)
                CF=N_factor
            if p_error<1e-15:
                p_error=0
            row.append(p_error)
        data.append(row)
    df = pd.DataFrame(data, index=CF, columns=diff_values)
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
    output.allow_fly = reactive.Value(False)
    output.allow_robot = reactive.Value(False)
    output.allow_decoder = reactive.Value(False)
    output.reply_green = reactive.Value(False)
    output.precomp_N=reactive.Value(0)
    output.precomp_diff=reactive.Value(0)
    # Track last submitted n_samp and differentiate
    output.last_submitted_n_samp = reactive.Value(0)
    output.last_submitted_differentiate = reactive.Value(0)
    output.database_reply_auto=reactive.Value("Waiting for pooling strategy file")
    output.biomek_csv=reactive.Value(0)
    output.hamilton_csv=reactive.Value(0)
    output.opentron_csv=reactive.Value(0)
    output.database_reply_decoder=reactive.Value("Waiting for pooling strategy file and readout")
    output.decoder_text=reactive.Value(0)
    output.decoder_diff=reactive.Value(-1)








    ls_met=['Pooling strategy', 'mean_experiments', 'max_compounds_per_well', 'n_wells', 'percentage_check', 'mean_extra_exp']
    output.summary_table = reactive.Value(pd.DataFrame(columns=ls_met))



    @reactive.Effect
    @reactive.event(input.submit)
    def _():
        # Get user inputs
        n_samp = input.n_samp()
        differentiate = input.differentiate()
        # Store last submitted values for use in outputs
        output.last_submitted_n_samp.set(n_samp)
        output.last_submitted_differentiate.set(differentiate)
        #print(n_samp)
        n_samp=output.last_submitted_n_samp.get()
        differentiate=output.last_submitted_differentiate.get()
        last_values_text = f"Max. number of samples: {n_samp}, Max. positives: {differentiate}"
        
        # Set output to display last submitted values
        output.last_values.set(last_values_text)

        if differentiate > n_samp:
            output_text=f'Maximum number of positives ({differentiate}) must always be smaller than the total number of samples ({n_samp})'
            output.database_reply.set(output_text)
            output.extra_computation.set(1)
            output.summary_table.set(pd.DataFrame(columns=ls_met))
            output.allow_downloads.set(False)
            output.allow_fly.set(False)
            output.reply_green.set(False)

        elif differentiate > MAX_DIFFERENTIATE:
            output_text_1=f'Maximum number of positives ({differentiate}) too high. The precomputed maximum is {MAX_DIFFERENTIATE}.\n'
            output_text_2 = '<span style="color: #228B22;">To download designs calculated on the fly for your specific values (or locally run the code), follow the section below.</span>'
            output_text=output_text_1+output_text_2
            output_text = output_text.replace('\n', '<br>')
            output.database_reply.set(output_text)
            output.extra_computation.set(1)
            output.summary_table.set(pd.DataFrame(columns=ls_met))
            output.allow_downloads.set(False)
            output.allow_fly.set(True)
            output.reply_green.set(False)

        elif n_samp > MAX_N:
            output_text_1=f'Maximum number of samples ({n_samp}) too high. The precomputed maximum is {MAX_N}.\n'
            output_text_2 = '<span style="color: #228B22;">To download designs calculated on the fly for your specific values (or locally run the code), follow the section below.</span>'
            output_text=output_text_1+output_text_2
            output_text = output_text.replace('\n', '<br>')
            output.database_reply.set(output_text)
            output.extra_computation.set(1)
            output.summary_table.set(pd.DataFrame(columns=ls_met))
            output.allow_downloads.set(False)
            output.allow_fly.set(True)
            output.reply_green.set(False)

        else:
            try:
                app_root = os.listdir('.')[0]
                #os.path.dirname(os.path.abspath(__file__))
                #os.listdir('.')[0]#next(f for f in os.listdir('.') if os.path.isdir(f) and f.startswith('precomp'))#next(f for f in os.listdir('.') if os.path.isdir(f) and not f.startswith('.'))#os.listdir('.')[0]
                WA_DIRECTORY = os.path.join(app_root, WA_SUB_DIRECTORY)
                result = find_closest_n_diff_folder(WA_DIRECTORY, n_samp, differentiate)
                if result:
                    folder_path, new_n, new_diff = result
                    #print(folder_path)
                    # Now folder_path is the path to the correct N_X/diff_Y folder
                    excel_filename = f'Metrics_N_{new_n}_diff_{new_diff}.csv'
                    excel_path = os.path.join(folder_path, excel_filename)
 
                    output.allow_downloads.set(True)
                    output.allow_fly.set(True)
                    if new_n!=n_samp or new_diff!=differentiate:
                        output_text_1=f'There is no precomputed summary table for {n_samp} samples with up to {differentiate} positives.\n'
                        output_text_2=f'The closest precomputed summary table is for {new_n} samples with up to {new_diff} positives, as shown below.\n'
                        output_text_3='<span style="color: #228B22;">To download designs calculated on the fly for your specific values (or locally run the code), follow the section below.</span>'
                        output_text=output_text_1+output_text_2+output_text_3
                        output_text = output_text.replace('\n', '<br>')
                        output.database_reply.set(output_text)
                        output.precomp_N.set(new_n)
                        output.precomp_diff.set(new_diff)
                        output.reply_green.set(False)
                        
                    else: 
                        output_text_1=f'There is a precomputed summary table for {n_samp} samples with up to {differentiate} positives.\n'
                        output_text_2 = '<span style="color: #228B22;">To download designs calculated on the fly for your specific values (or locally run the code), follow the section below.</span>'
                        output_text=output_text_1+output_text_2
                        output_text = output_text.replace('\n', '<br>')
                        output.database_reply.set(output_text)
                        output.precomp_N.set(n_samp)
                        output.precomp_diff.set(differentiate)
                        output.reply_green.set(True)



                        #print(output_text)
                    if os.path.isfile(excel_path):
                        metrics_data = pd.read_csv(excel_path)

                        # Use metrics_data as needed
                    else:
                        # Handle missing metrics file
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
                dict_ren={'N wells':'N pools (W)', 'Max compunds per well': 'Max compounds in a pool'}#{i:j for i,j in zip(metrics_data.columns,ls_met_nice)}
                metrics_data.rename(columns=dict_ren, inplace=True)
                
                #metrics_data.rename(index={'Chinese trick':'Chinese reminder'}, inplace=True)
                output.summary_table.set(metrics_data)

                copio=metrics_data.copy()
                if 'Mean experiments' in copio.columns:
                    copio = copio[copio['Mean experiments'] < n_samp]

                #long_exps=f'If your experiments are <b>expensive</b>, the best method for you might be the <b>{copio.iloc[0,0]}</b>.<br>'
                if 'Mean experiments' in copio.columns and 'Mean steps' in copio.columns:
                    copio = copio.sort_values(by=['Mean experiments', 'Mean steps'], ascending=[True, True])
                    long_exps=f'If you mainly need to <b>reduce test number</b>, the best method for you might be the <b>{copio.iloc[0,0]}</b>.<br>'

                #expensive_exps=f'If your experiments are <b>long</b>, the best method for you might be the <b>{copio.iloc[0,0]}</b>.<br>'
                if 'Mean experiments' in copio.columns and 'Mean steps' in copio.columns:
                    copio = copio.sort_values(by=['Mean steps', 'Mean experiments'], ascending=[True, True])
                    expensive_exps=f'If you need results <b>quickly</b>, the best method for you might be the <b>{copio.iloc[0,0]}</b>.<br>'
                
                #low_pools=f'If you <b>cannot pool many samples</b> together, the best method for you might be the <b>{copio.iloc[0,0]}</b>.<br>'
                if 'Mean experiments' in copio.columns and 'Max samples per pool' in copio.columns:
                    copio = copio.sort_values(by=['Max samples per pool', 'Mean experiments'], ascending=[True, True])                
                    low_pools=f'If you <b>cannot pool many samples</b> together, the best method for you might be the <b>{copio.iloc[0,0]}</b>.<br>'
                

                
                #long_exps_html = "<br>".join([f"<span style='font-size:15px; color:#444;'>{col}</span>" for col in copio.columns.tolist()])

                output.summary_text.set(expensive_exps+long_exps+low_pools)#+long_exps_html)

                #TBLS=CR[1]
                DFFS={}
                table_path=os.path.join(folder_path,'WAs')
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
                app_root =os.listdir('.')[0]# next(f for f in os.listdir('.') if os.path.isdir(f) and f.startswith('precomp'))#
                print(files)
                print(os.getcwd())
                output.debug.set(str([f for f in files]))
                output.debug.set(f"cwd: {os.getcwd()}, files: {os.listdir('.')}")

                app_root =os.listdir('.')[0]# next(f for f in os.listdir('.') if os.path.isdir(f) and not f.startswith('.'))#os.listdir('.')[0]  # first (and usually only) item — the app folder
                app_root_contents = os.listdir(app_root)

                if 'precomputed' in app_root_contents:
                    output.debug.set(f"precomputed/ FOUND in {app_root}")
                    precomputed_path = os.path.join(app_root, 'precomputed')
                    output.debug.set(f"Files in precomputed/: {os.listdir(precomputed_path)}")
                else:
                    output.debug.set(f"precomputed/ NOT found in {app_root}. Contents: {app_root_contents}")
                #diff_folder, new_diff = find_closest_diff_folder(n_folder_path, differentiate)
                output_text_1 = '<span style="color: #c00;">Something went wrong while looking for your file.</span> <br>'
                output_text_2 = '<span style="color: #228B22;">To download designs calculated on the fly for your specific values (or locally run the code), follow the section below.</span>'
                output_text=output_text_1+output_text_2
                output_text = output_text.replace('\n', '<br>')
                output.database_reply.set(output_text)
                output.summary_table.set(pd.DataFrame(columns=ls_met))
                output.allow_downloads.set(False)
                
            
            
        # Prepare the correct command for pool_N.py based on its arguments
        command_p = f"python pool_N.py --n_samp {n_samp} --differentiate {differentiate} --method all --path your/path"
        intro_command='You can also run the code locally to generate pooling designs. For the values selected above, run:\n'
        mid_github_command='(see GitHub for details)\n\n'
        output.personalized_command.set(intro_command+mid_github_command+command_p)


        intro_command='To decode the results of your tests done with a designs from PoolPy, run:\n'
        mid_github_command='(see GitHub for details)\n\n'
        decode_p = f"python decode_N.py --differentiate {differentiate} path_to_WA your/path/to/pooling_strategy.csv --readout path/to/readout.csv (or your readout)\n\n"
        readout_exp='Your readout should either be in binary form (e.g . 1,0,0,1,0) or as a list of positive pools (e.g. 0,3).'
        output.decoder.set(intro_command+decode_p+mid_github_command+readout_exp)

            #md=np.max(np.array(list(f2)))

            


    @reactive.Effect
    @reactive.event(input.prev_submit)
    def _():
        n = input.prev_n_samp()
        p = input.prev_prevalence()
        max_e = input.prev_max_error()
        output.prev_n_samp.set(n)
        output.prev_prevalence.set(p)
        output.prev_max_error.set(max_e)
        n = output.prev_n_samp.get()
        p = output.prev_prevalence.get()
        max_e = output.prev_max_error.get()
        output.last_val_prev.set(f"Number of Samples: {n}, Prevalence: {p}, Max. error: {max_e}")

    @reactive.Effect
    @reactive.event(input.uploaded_csv_auto)
    def _():
        fileinfo = input.uploaded_csv_auto()
        if fileinfo is not None:
            if fileinfo[0]["size"] > 100 * 1024:
                output.database_reply_auto.set("File too large. Please upload a CSV smaller than 100kb.")
                output.allow_robot.set(False)
            else:
                output.database_reply_auto.set("File is ok, computing automated pooling scripts.")
                output.allow_robot.set(True)
                WA_df = parsed_file_auto()
                WA=WA_df.values
                offsetted=int((WA.shape[0]//96)+1)
                ls_to_df=[]
                for i,j in zip(*np.where(WA==1)):
                    d_plate, d_well= from_id_to_well(i)
                    r_plate, r_well= from_id_to_well(j, offset=offsetted)
                    ls_to_df.append([d_plate,d_well,r_plate,r_well])

                general_df_out=pd.DataFrame(data=ls_to_df, columns=['Source', 'SourceWell', 'Dest', 'DestWell'])
                general_df_out['Volume']=1

                biomek_out=general_df_out.copy()
                biomek_out['Source'] = biomek_out['Source'].apply(lambda x: f'Plate{x}')
                biomek_out['Dest'] = biomek_out['Dest'].apply(lambda x: f'Plate{x}')
                output.biomek_csv.set(biomek_out)

                hamilton_out=general_df_out.copy()
                hamilton_out.rename(columns={'Source':'SourceSite', 'Dest':'TargetSite', 'DestWell':'TargetWell'}, inplace=True)
                output.hamilton_csv.set(hamilton_out)

                opentron_out=general_df_out.copy()
                opentron_out.rename(columns={'Source':'Source Slot','SourceWell':'Source Well', 'Dest':'Dest Slot', 'DestWell':'Dest Well', 'Volume': 'Volume (in ul)'}, inplace=True)
                opentron_out['Source Labware']=opentron_out['Source Slot'].apply(lambda x : f'96_wells_plate_{x}')
                opentron_out['Dest Labware']=opentron_out['Dest Slot'].apply(lambda x : f'96_wells_plate_{x}')
                opentron_out['Source Aspiration Height Above Bottom (in mm)']=1
                opentron_out=opentron_out[['Source Labware', 'Source Slot', 'Source Well', 'Source Aspiration Height Above Bottom (in mm)', 'Dest Labware', 'Dest Slot', 'Dest Well', 'Volume (in ul)']]
                output.opentron_csv.set(opentron_out)


                output.database_reply_auto.set("Automated pooling scripts ready.")

                
                

    @reactive.calc
    def parsed_file_auto():
        file = input.uploaded_csv_auto()
        if file is None:
            return pd.DataFrame()
        return pd.read_csv(file[0]["datapath"], index_col=0)
    
    @reactive.calc
    def parsed_file_decoder():
        file = input.uploaded_csv_decoder()
        if file is None:
            return pd.DataFrame()
        return pd.read_csv(file[0]["datapath"], index_col=0)

    @reactive.Effect
    @reactive.event(input.readout_submit)
    def _():
        fileinfo = input.uploaded_csv_decoder()
        output.database_reply_decoder.set("Please upload a pooling strategy file.")
        if fileinfo is not None:
            if fileinfo[0]["size"] > 100 * 1024:
                output.database_reply_decoder.set("File too large. Please upload a CSV smaller than 100kb.")
                output.allow_decoder.set(False)
            elif fileinfo[0]["size"] < 100 * 1024:
                file_name = fileinfo[0]["name"] if fileinfo and "name" in fileinfo[0] else "No file uploaded"
                #output.database_reply_decoder.set("File and readout are ok, decoding.")
                output.allow_decoder.set(True)
                readout_str = input.readout_string()
                # Now readout_str contains the string entered by the user, e.g. "0,1,1,0,0,0,1" or "3,7,4"
                # You can process it as needed:
                try:
                    readout_list = []
                    for x in readout_str.split(','):
                        x_clean = x.strip()
                        if x_clean == '':
                            continue
                        if not x_clean.isdigit():
                            msg = f"<b>Error:</b> Invalid entry '{x_clean}' in readout. Please enter only integers separated by commas.<br>"
                            output.database_reply_decoder.set(msg)
                            return
                        readout_list.append(int(x_clean))
                except Exception as e:
                    msg = f"<b>Error:</b> Could not parse readout. Please enter only integers separated by commas.<br>"
                    output.database_reply_decoder.set(msg)
                    return
                diff_deco = input.decoder_diff()
                try:
                    diff_deco = int(diff_deco)
                except :
                    diff_deco = -1

                output.decoder_diff.set(diff_deco)
                diff_deco=output.decoder_diff.get()

                msg=f'Processing file <b>{fileinfo[0]["name"]}</b> <br>'

                if len(readout_list)>0:
                    WA_df = parsed_file_decoder()
                    WA=WA_df.values
                    n_pools=WA.shape[1]
                    n_compounds=WA.shape[0]
                    if diff_deco==-1:
                        diff_deco=n_compounds

                    readout_list.sort()
                    msg=f'Processing file <b>{fileinfo[0]["name"]}</b> with max. <b>{diff_deco}</b> positive samples and readout: <br>'
                    msg+=f"{readout_list}<br>"
                    #msg+=f'{WA}'
                    readout=np.array(readout_list, dtype=int)
                    if max(readout)>n_pools:
                        invalid_entries = [val for val in readout_list if val > n_pools]
                        msg += f"<br>The following entries in your readout are bigger than the number of pools ({n_pools}):<br> <b>{invalid_entries}</b>"
                        msg+="<br>Please <b>correct your readout</b>.<br>"

                        output.database_reply_decoder.set(msg)
                        return
                    if np.max(readout)>1 or len(readout)!=n_pools:
                        readout_bin_ls = [1 if i in readout else 0 for i in range(n_pools)]
                        readout=np.array(readout_bin_ls)
                    msg+=f'<br>The uploaded pooling strategy comprizes {n_compounds} samples in {n_pools} pools.<br><br>'
                    readout_bl=np.array(readout.astype(bool).astype(int))
                    mask = ~np.any((WA == 1) & (readout_bl == 0), axis=1)
                    #msg+=f'boolean readout {readout_bl}<br>'
                    original_indices = np.where(mask)[0]  # Get original row indices
                    filtered_WA = WA[mask]
                    n_compounds=filtered_WA.shape[0]

                    if diff_deco> n_compounds:
                        diff_deco=n_compounds
                    #msg+=f'<br>we have {filtered_WA.shape} shape of the filtered WA<br>'

                    if n_compounds<2:
                        if n_compounds==1:
                            decoded=[original_indices[0]]
                            msg += 'A single positive sample was found:<br>'
                            for deco in decoded:
                                msg += f'<b>Sample: {deco}</b><br>'
                        else:
                            msg += '<b>We found no matches for the given parameters, check your input or try increasing the differentiate value.</b>'
                    
                    else:
                        scrambler={1:np.arange(n_compounds)}
                        for j in range(2,diff_deco+1):
                            scrambler.update({j:np.array(list(itertools.combinations(np.arange(n_compounds),j)))})
                        decoded_pre=decode_precomp(well_assigner=filtered_WA, differentiate= diff_deco, scrambler=scrambler, 
                                readout=readout_bl)
                        
                        # Map filtered indices back to original indices
                        decoded = [[int(original_indices[idx]) for idx in combination] for combination in decoded_pre]
                        
                        # Remove duplicate combinations (convert to tuples for uniqueness, then back to lists)

                        if len(decoded)==0:
                            msg+= '<b>We found no matches for the given parameters, check your input or try increasing the differentiate value.'

                        elif len(decoded)==1:
                            msg+=f'A single possible combination of positive samples was found. The positive samples are:<br>'
                            msg += f'<b>Samples: {decoded[0]}</b><br>'

                        elif len(decoded)>n_compounds:
                            decoded_set = list(set([x for combo in decoded for x in combo]))
                            decoded_set.sort()
                            msg += (
                                '<b>Putative positive samples were identified</b>, but the exact combination could not be pinpointed.<br>'
                                'Either test all putative positive samples individually or change pooling strategy. A lower differentiate (only if it makes sense) might narrow it down.<br>'
                                f'There are up to {min([diff_deco,len(decoded_set)])} positive samples among the following samples: <b>{decoded_set}</b>.'
                            )
                        else:
                            msg += f'{len(decoded)} possible combinations of positive samples were found. The possible combinations are:<br>'
                            for i,deco in enumerate(decoded):
                                if i!=0:
                                    msg+='or<br>'
                                msg += f'Samples <b>{deco}</b><br>'

                    output.database_reply_decoder.set(msg)

                    mess = output.database_reply_decoder.get()
                    message = re.sub(r'</?b>', '', mess)
                    txt = f"Uploaded file: {file_name}\nReadout: {readout_str}\n\n{message.replace('<br>', '\n')}"
                    output.decoder_text.set(txt)

                else:
                    msg+='Please provide a readout in the correct format.'
                    output.database_reply_decoder.set(msg)


            #elif False:
           


    @output
    @render.ui
    def last_val_prev():
        val = output.last_val_prev.get()
        if not val:
            return ui.HTML("")
        html = f'<div style="font-weight: bold; text-align: center; margin: 18px 0 18px 0; font-size: 17px;">{val}</div>'
        return ui.HTML(html)


    @output
    @render.text
    def allow_downloads():
        return str(output.allow_downloads.get()).lower()  # returns 'true' or 'false'
    
    @output
    @render.text
    def allow_fly():
        return str(output.allow_fly.get()).lower() 
    @output
    @render.text
    def allow_robot():
        return str(output.allow_robot.get()).lower() 
    
    @output
    @render.text
    def allow_decoder():
        return str(output.allow_decoder.get()).lower() 

    @output
    @render.ui
    def last_val():
        val = output.last_values.get()
        if not val:
            return ui.HTML("")
        html = f'<div style="font-weight: bold; text-align: center; margin: 18px 0 18px 0; font-size: 17px;">{val}</div>'
        return ui.HTML(html)

    @output
    @render.ui
    def database_r():
        reply = output.database_reply.get()
        # Determine color: green if reply contains 'precomputed strategy', else red
        color = '#2ecc40' if output.reply_green.get() else '#ff4136'
        html = f'<div style="text-align: center; margin-top: 18px;">'
        if reply:
            html += f'<div style="color: {color}; font-size: 20px; font-weight: bold;">{reply}</div>'
        html += '</div>'
        return ui.HTML(html)
    
    @output
    @render.ui
    def database_r_auto():
        reply = output.database_reply_auto.get()
        html = f'<div style="text-align: center; margin-top: 18px;">'
        if reply:
            html += f'<div style="color: #111; font-size: 20px; font-weight: normal;">{reply}</div>'
        html += '</div>'
        return ui.HTML(html)
    
    @output
    @render.ui
    def database_r_decoder():
        reply = output.database_reply_decoder.get()
        html = f'<div style="text-align: center; margin-top: 18px;">'
        if reply:
            html += f'<div style="color: #111; font-size: 20px; font-weight: normal;">{reply}</div>'
        html += '</div>'
        return ui.HTML(html)
    
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


# On the fly downloads:
    @output
    @render.download(filename=lambda: "matrix_pooling.csv")
    async def download_table_matrix_fly():
        # Yield the content of the CSV file
        DFF=assign_wells_mat(n_compounds=output.last_submitted_n_samp.get())
        DFFS=clean_WA(DFF)
        yield DFFS.to_csv(index=True)
    
    @output
    @render.download(filename=lambda: "2D_pooling.csv")
    async def download_table_2d_fly():
        # Yield the content of the CSV file
        DFF=assign_wells_multidim(n_compounds=output.last_submitted_n_samp.get(), n_dims=2)
        DFFS=clean_WA(DFF)
        yield DFFS.to_csv(index=True)

    @output
    @render.download(filename=lambda: "3D_pooling.csv")
    async def download_table_3d_fly():
        # Yield the content of the CSV file
        DFF=assign_wells_multidim(n_compounds=output.last_submitted_n_samp.get(), n_dims=3)
        DFFS=clean_WA(DFF)
        yield DFFS.to_csv(index=True)

    @output
    @render.download(filename=lambda: "4D_pooling.csv")
    async def download_table_4d_fly():
        # Yield the content of the CSV file
        DFF=assign_wells_multidim(n_compounds=output.last_submitted_n_samp.get(), n_dims=4)
        DFFS=clean_WA(DFF)
        yield DFFS.to_csv(index=True)

    @output
    @render.download(filename=lambda: "5D_pooling.csv")
    async def download_table_5d_fly():
        # Yield the content of the CSV file
        DFF=assign_wells_multidim(n_compounds=output.last_submitted_n_samp.get(), n_dims=5)
        DFFS=clean_WA(DFF)
        yield DFFS.to_csv(index=True)

    @output
    @render.download(filename=lambda: "6D_pooling.csv")
    async def download_table_6d_fly():
        # Yield the content of the CSV file
        DFF=assign_wells_multidim(n_compounds=output.last_submitted_n_samp.get(), n_dims=6)
        DFFS=clean_WA(DFF)
        yield DFFS.to_csv(index=True)

    @output
    @render.download(filename=lambda: "binary_pooling.csv")
    async def download_table_binary_fly():
        # Yield the content of the CSV file
        DFF=assign_wells_bin(n_compounds=output.last_submitted_n_samp.get())
        DFFS=clean_WA(DFF)
        yield DFFS.to_csv(index=True)

    @output
    @render.download(filename=lambda: "STD_pooling.csv")
    async def download_table_STD_fly():
        # Yield the content of the CSV file
        DFF=assign_wells_STD(n_compounds=output.last_submitted_n_samp.get(), differentiate=output.last_submitted_differentiate.get())
        DFFS=clean_WA(DFF)
        yield DFFS.to_csv(index=True)

    @output
    @render.download(filename=lambda: "Chinese_remainder_pooling.csv")
    async def download_table_CT_fly():
        # Yield the content of the CSV file
        DFF=assign_wells_chinese(n_compounds=output.last_submitted_n_samp.get(), differentiate=output.last_submitted_differentiate.get())
        DFFS=clean_WA(DFF)
        yield DFFS.to_csv(index=True)

    @output
    @render.download(filename=lambda: "Chinese_remainder_bktrk_pooling.csv")
    async def download_table_CT_bktrk_fly():
        # Yield the content of the CSV file
        DFF=assign_wells_chinese(n_compounds=output.last_submitted_n_samp.get(), differentiate=output.last_submitted_differentiate.get(), backtrack=True)
        DFFS=clean_WA(DFF)
        yield DFFS.to_csv(index=True)

    @output
    @render.download(filename=lambda: "Chinese_remainder_special_pooling.csv")
    async def download_table_CT_special_fly():
        # Yield the content of the CSV file
        if output.last_submitted_differentiate.get()==2 or output.last_submitted_differentiate.get()==3:
            DFF=assign_wells_chinese(n_compounds=output.last_submitted_n_samp.get(), differentiate=output.last_submitted_differentiate.get(), special_diff=True)
            DFFS=clean_WA(DFF)
            yield DFFS.to_csv(index=True)
        else:
            DFF=assign_wells_chinese(n_compounds=output.last_submitted_n_samp.get(), differentiate=output.last_submitted_differentiate.get(), backtrack=True)
            DFFS=clean_WA(DFF)
            yield DFFS.to_csv(index=True)

    @output
    @render.ui
    def fly_download_ind_title():
        n_samp = output.last_submitted_n_samp.get()
        differentiate = output.last_submitted_differentiate.get()
        return ui.HTML(
            f'<div style="margin-bottom: 20px; font-size: 22px;">'
            f'<b>Differentiate-indipendent</b> designs calculated on the fly for {n_samp} samples with up to {differentiate} positives'
            f'<br> <b style="color: #c00;">Use with caution for differentiate &gt;1</b>'
            f'</div>'
        )
        
    @output
    @render.ui
    def fly_download_dip_title():
        n_samp = output.last_submitted_n_samp.get()
        differentiate = output.last_submitted_differentiate.get()
        return ui.HTML(
            f'<div style="margin-bottom: 20px; font-size: 22px;">'
            f'<b>Differentiate-dependent</b> designs calculated on the fly for {n_samp} samples with up to {differentiate} positives'
            f'</div>'
        )


    @output
    @render.image
    def logo_img_main():
        here = Path(__file__).parent
        return {
            "src": here / "static/banner_noborder.png",
            "style": (
                "height: 130px; "
                "width: 100vw; "
                "object-fit: cover; "
                "display: block; "
                "margin-left: 50%; "
                "transform: translateX(-50%); "
                "overflow: hidden; "
                "padding: 0;"
            )
        }
    
    @output
    @render.image
    def logo_img_methods():
        here = Path(__file__).parent
        return {
            "src": here / "static/banner_noborder.png",
            "style": (
                "height: 130px; "
                "width: 100vw; "
                "object-fit: cover; "
                "display: block; "
                "margin-left: 50%; "
                "transform: translateX(-50%); "
                "overflow: hidden; "
                "padding: 0;"
            )
        }

    @output
    @render.image
    def logo_img_prev():
        here = Path(__file__).parent
        return {
            "src": here / "static/banner_noborder.png",
            "style": (
                "height: 130px; "
                "width: 100vw; "
                "object-fit: cover; "
                "display: block; "
                "margin-left: 50%; "
                "transform: translateX(-50%); "
                "overflow: hidden; "
                "padding: 0;"
            )
        }
    
    @output
    @render.image
    def logo_img_automation():
        here = Path(__file__).parent
        return {
            "src": here / "static/banner_noborder.png",
            "style": (
                "height: 130px; "
                "width: 100vw; "
                "object-fit: cover; "
                "display: block; "
                "margin-left: 50%; "
                "transform: translateX(-50%); "
                "overflow: hidden; "
                "padding: 0;"
            )
        }

    @output
    @render.image
    def logo_img_decoder():
        here = Path(__file__).parent
        return {
            "src": here / "static/banner_noborder.png",
            "style": (
                "height: 130px; "
                "width: 100vw; "
                "object-fit: cover; "
                "display: block; "
                "margin-left: 50%; "
                "transform: translateX(-50%); "
                "overflow: hidden; "
                "padding: 0;"
            )
        }

    @output
    @render.ui
    def summary_table_title():
        n_samp = int(output.last_submitted_n_samp.get())
        diff = int(output.last_submitted_differentiate.get())
        new_n = output.precomp_N.get()
        new_diff = output.precomp_diff.get()
        allow_downloads = output.allow_downloads.get()
        # Only show the 'closest precomputed values' message if user has submitted and allow_downloads is True
        if allow_downloads:
            if new_diff != diff or new_n != n_samp:
                return ui.HTML(f'<div style="text-align: center; font-weight: bold; color: #c00; font-size: 22px;">Summary pooling table for closest precomputed values (S={new_n}, D={new_diff})</div>')
            else:
                return ui.h4("Summary pooling strategy table", style="text-align: center;")
        else:
            return ui.h4("Summary pooling strategy table", style="text-align: center;")

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
        n_samp = int(output.last_submitted_n_samp.get())
        # Color logic
        green_indices = {}
        # Green for three smallest in 'Mean experiments' and 'Max compounds in a pool',
        # but only if value is less than n_samp (do not color green if value >= n_samp)
        for col in ['Mean experiments', 'Max compounds in a pool']:
            if col in df.columns:
                col_vals = df[col]
                idxs = pd.Series(col_vals).nsmallest(3).index.tolist()
                # Remove any index where value >= n_samp
                idxs = [i for i in idxs if col_vals[i] < n_samp]
                green_indices[col] = idxs

        # Remove green logic for 'N pools' entirely (do not color green in N pools)
        # (No code here)

        # Red or green for 'Mean steps' (if >3, color all ==1 green)
        red_indices = []
        if 'Mean steps' in df.columns:
            col_vals = df['Mean steps']
            gt1 = col_vals[col_vals > 1]
            if 1 <= len(gt1) <= 3:
                red_indices = gt1.index.tolist()
            elif len(gt1) > 3:
                green_indices['Mean steps'] = col_vals.index[col_vals == 1].tolist()

        for i, row in df.iterrows():
            # Robustly get the mean experiments value for this row
            row_red = False
            mean_exp = None
            mean_exp_col = None
            for c in df.columns:
                if c.strip().lower().replace(' ', '') == 'meanexperiments':
                    mean_exp_col = c
                    break
            if mean_exp_col and n_samp is not None:
                try:
                    mean_exp = float(row[mean_exp_col])
                    if mean_exp > n_samp:
                        row_red = True
                except Exception:
                    pass
            html += '<tr>'
            for col in df.columns:
                val = row[col]
                td_style = 'border: 1px solid #aaa; padding: 8px 12px; text-align: center;'
                if col == 'N pools':
                    td_style += ' min-width: 120px;'
                # Red for full row if Mean experiments > n_samp
                if row_red:
                    td_style += ' background-color: #ffcccc; color: #111;'
                # Green logic (but NOT for method if row_red)
                if col in green_indices and i in green_indices[col]:
                    if not row_red:
                        td_style += ' background-color: #b6fcb6; color: #111;'
                # Red logic for Mean steps (overrides green)
                if col == 'Mean steps' and i in red_indices:
                    td_style += ' background-color: #ffcccc; color: #111;'
                # Debug: show mean_exp value in a tooltip for the method column
                if col.strip().lower().replace(' ', '') == 'poolingstrategy':
                    html += f'<td style="{td_style}" title="mean_exp={mean_exp}">{val}</td>'
                else:
                    html += f'<td style="{td_style}">{val}</td>'
            html += '</tr>'
        html += '</table>'
        return ui.HTML(html)
  

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
        html += '<tr><th style="border: 1px solid #aaa; padding: 10px 16px; text-align: center;">S \\ D</th>'
        for col in df_orig.columns:
            html += f'<th style="border: 1px solid #aaa; padding: 10px 16px; text-align: center;">{col}</th>'
        html += '</tr>'
        for i, (idx, row) in enumerate(df_orig.iterrows()):
            html += f'<tr><th style="border: 1px solid #aaa; padding: 10px 16px; text-align: center;">{idx}</th>'
            for j, val in enumerate(row):
                # Compare by position, not by index/column label
                try:
                    val_corr = df_corr.iloc[i, j]
                except Exception:
                    val_corr = None
                style = color_for_value(val, val_corr)
                html += f'<td style="border: 1px solid #aaa; padding: 10px 16px; text-align: center; {style}">{val:.3g}</td>'
            html += '</tr>'
        html += '</table>'
        html += '</div>'

        # Corrected table
        html += '<div>'
        html += '<div style="text-align: center; font-weight: bold; margin-bottom: 8px;">Corrected for FWER</div>'
        html += '<table style="border-collapse: collapse; font-size: 18px; min-width: 300px;">'
        html += '<tr><th style="border: 1px solid #aaa; padding: 10px 16px; text-align: center;">S \\ D</th>'
        for col in df_corr.columns:
            html += f'<th style="border: 1px solid #aaa; padding: 10px 16px; text-align: center;">{col}</th>'
        html += '</tr>'
        for i, (idx, row) in enumerate(df_corr.iterrows()):
            html += f'<tr><th style="border: 1px solid #aaa; padding: 10px 16px; text-align: center;">{idx}</th>'
            for j, val_corr in enumerate(row):
                try:
                    val_orig = df_orig.iloc[i, j]
                except Exception:
                    val_orig = None
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
        html += '<div style="margin-bottom: 6px;"><span style="background-color: #ffeaea; color: #111; padding: 4px 14px; border-radius: 4px; border: 1px solid #eee;">Red: pooling error &gt; max. error for both single pooling and FWER (all poolings combined).</span></div>'
        html += '<div style="margin-bottom: 6px;"><span style="background-color: #ffffea; color: #111; padding: 4px 14px; border-radius: 4px; border: 1px solid #eee;">Yellow: pooling error ≤ max. error for single pooling, but &gt; max. error for FWER (all poolings combined).</span></div>'
        html += '<div><span style="background-color: #b6fcb6; color: #111; padding: 4px 14px; border-radius: 4px; border: 1px solid #eee;">Green: pooling error ≤ max. error for both single pooling and FWER (all poolings combined).</span></div>'
        html += '</div>'
        return ui.HTML(html)

    @output
    @render.ui
    def prevalence_explanation():
        html = '<div style="margin-top: 28px; font-size: 14px; color: #444; text-align: center;">These tables are meant to help decide on a pooling strategy for specific use cases.</div>'
        return ui.HTML(html)
    

    @output
    @render.download(filename=lambda: "biomek_robot_code.csv")
    async def download_biomek():
        df_out=output.biomek_csv.get()
        yield df_out.to_csv(index=False)

    @output
    @render.download(filename=lambda: "hamilton_robot_code.csv")
    async def download_hamilton():
        df_out=output.hamilton_csv.get()
        yield df_out.to_csv(index=False)

    @output
    @render.download(filename=lambda: "opentron_robot_code.csv")
    async def download_opentron():
        df_out=output.opentron_csv.get()
        yield df_out.to_csv(index=False)



    @output
    @render.download(filename=lambda: "decoder_output.txt")
    async def download_decoder_output():
        text=output.decoder_text.get()

        yield text

app = App(app_ui, server)


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
