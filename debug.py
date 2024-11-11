from shiny import App, reactive, render, ui

# Define a simple function to calculate the square of a number
def square(value):
    return value ** 2

# Define the UI
app_ui = ui.page_fluid(
    # Input field to enter a number
    ui.input_numeric("input_value", "Enter a number:", value=1),
    
    # Button to trigger the square calculation
    ui.input_action_button("submit", "Calculate Square"),
    
    # Output area to display the result
    ui.h4("Result:"),
    ui.output_text_verbatim("output_square")  # Placeholder for output
)

# Define the server logic
def server(input, output, session):

    # Reactive value to store the result
    result = reactive.Value(0)

    @reactive.Effect
    @reactive.event(input.submit)  # Trigger only when 'Submit' button is clicked
    def _():
        # Get input value from user and calculate its square
        entered_value = input.input_value()
        squared_value = square(entered_value)

        # Print for debugging (optional)
        print(f"Entered Value: {entered_value}, Square: {squared_value}")

        # Set result to display in UI
        result.set(f"The square of {entered_value} is {squared_value}.")

    # Define a reactive text output using @render.text
    @output
    @render.text
    def output_square():
        return result.get()

# Create the app object
app = App(app_ui, server)
