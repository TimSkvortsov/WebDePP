async function main() {
    disableAllButtons(true);
    setLoadingStatus('Loading packages and initialising the program, please wait...');
    document.getElementById('predict-depolymerases').disabled = true;
    let pyodide = await loadPyodide();
    await pyodide.loadPackage("micropip");
    const micropip = pyodide.pyimport("micropip");
    await micropip.install('https://timskvortsov.github.io/pyodide/depp-1.0.0-py3-none-any.whl');
    pyodide.runPython(`
        from depp import depp
        import js
    `);



    // Initialize global variables to store protein parameters and depolymerase predictions
    let protein_parameters;
    let depolymerase_predictions;

    // disableAllButtons(false);
    document.getElementById('reset').disabled = false;
    setLoadingStatus('Initialisation complete!', true);
    updateButtonStates();

    // Add event listeners for buttons and file inputs
    document.getElementById('fasta-upload').addEventListener('input', updateCalculateParametersButtonState);
    document.getElementById('calculate-parameters').addEventListener('click', calculateProteinParameters);
    document.getElementById('predict-depolymerases').addEventListener('click', predictPhageDepolymerases);
    document.getElementById('reset').addEventListener('click', reset);
    
    document.getElementById('display-predictions').addEventListener('click', () => {
        if (!depolymerase_predictions) {
            alert("Phage depolymerase predictions are not available. Please calculate them first.");
            return;
        }
        displayPredictionsTable();
    });
    
    document.getElementById('save-parameters').addEventListener('click', () => {
        if (!protein_parameters) {
            alert("Protein parameters are not available. Please calculate them first.");
            return;
        }
        saveJSONToCSV(protein_parameters, 'protein_parameters.csv');
    });
    
    document.getElementById('save-predictions').addEventListener('click', () => {
        if (!depolymerase_predictions) {
            alert("Phage depolymerases predictions are not available. Please calculate them first.");
            return;
        }
        saveJSONToCSV(depolymerase_predictions, 'depolymerase_predictions.csv');
    });
    
    
    async function calculateProteinParameters() {
        const calculateParametersButton = document.getElementById('calculate-parameters');
        toggleButtonLoading(calculateParametersButton, true);
        const fastaFileInput = document.getElementById('fasta-upload');
        const file = fastaFileInput.files[0];
        if (!file) {
            toggleButtonLoading(calculateParametersButton, false);
            alert('Please select a FASTA file.');
            return;
        }
        const fileContent = await file.text();
        try {
            // Set the file content to a Python variable
            pyodide.globals.set('fileContent', fileContent);
    
            // Run the Python code
            const protein_parameters_json = await pyodide.runPythonAsync(`
                import json
                from Bio import SeqIO
                from Bio.SeqIO.FastaIO import SimpleFastaParser
                from io import StringIO
                from Bio.SeqRecord import SeqRecord
                from Bio.Seq import Seq

                # Create a StringIO object with the file content
                fasta_file = StringIO(fileContent)

                # Read the FASTA records using the depp.fasta_read function
                seq_records = depp.read_fasta(fasta_file)

                # Calculate protein parameters
                protein_parameters = [depp.calculate_protein_parameters(seq_record) for seq_record in seq_records]
                protein_df = depp.generate_dataframe(protein_parameters)

                protein_df.to_json(orient='records')
            `);

            // Parse the JSON string and display the result in the console or process it further
            protein_parameters = JSON.parse(protein_parameters_json);
            toggleButtonLoading(calculateParametersButton, false);
            updateButtonStates();
            

        } catch (error) {
            console.log(error.message);
        } finally {
            toggleButtonLoading(calculateParametersButton, false);
            updateButtonStates();
        }
        
    }
    
    
    async function predictPhageDepolymerases() {
        const predictDepolymerasesButton = document.getElementById('predict-depolymerases');
        toggleButtonLoading(predictDepolymerasesButton, true);
        if (!protein_parameters) {
            alert("Protein parameters are not available. Please calculate them first.");
            toggleButtonLoading(predictDepolymerasesButton, false);
            return;
        }

        const trainingFileInput = document.getElementById('training-upload');
        const file = trainingFileInput.files[0];
        let trainingFileContent;

        if (file) {
            trainingFileContent = await file.text();
        } else {
            const response = await fetch('./TrainingSet/TrainingSet.csv');
            trainingFileContent = await response.text();
        }

        try {
            // Set the training file content and protein parameters to Python variables
            pyodide.globals.set('trainingFileContent', trainingFileContent);
            pyodide.globals.set('protein_parameters_json', JSON.stringify(protein_parameters));

            // Run the Python code to train the Random Forest model and make predictions
            const depolymerase_predictions_json = await pyodide.runPythonAsync(`
                import json
                import os
                import pandas as pd
                import numpy as np
                from sklearn.ensemble import RandomForestClassifier
                from sklearn.pipeline import Pipeline
                from sklearn.preprocessing import PolynomialFeatures, MinMaxScaler
                from io import StringIO

                # Save the training dataset to a temporary file
                temp_training_set_path = "temp_training_set.csv"
                with open(temp_training_set_path, "w") as temp_training_set_file:
                    temp_training_set_file.write(trainingFileContent)

                # Train the Random Forest model
                rf_model = depp.train_model(temp_training_set_path)

                # Remove the temporary file
                os.remove(temp_training_set_path)

                # Create a DataFrame from JSON
                protein_parameters_dict = json.loads(protein_parameters_json)
                protein_df = pd.DataFrame.from_dict(protein_parameters_dict)

                # Predict phage depolymerases
                depolymerase_predictions = depp.predict_depolymerases(protein_df, rf_model)

                # Serialize the predictions to JSON
                depolymerase_predictions.to_json(orient='records')
            `);

            // Parse the JSON string and display the result in the console or process it further
            depolymerase_predictions = JSON.parse(depolymerase_predictions_json);
            toggleButtonLoading(predictDepolymerasesButton, false);
            updateButtonStates();

        } catch (error) {
            console.log(error.message);
        } finally {
            toggleButtonLoading(predictDepolymerasesButton, false);
            updateButtonStates();
        }
    }


    function JSONtoCSV(jsonData) {
        // Convert JSON data to an array of objects
        const data = Object.keys(jsonData).map(key => jsonData[key]);
    
        // Extract the header (keys) from the first object
        const header = Object.keys(data[0]);
    
        // Convert JSON data to CSV
        const csvContent =
            header.join(',') + '\n' +
            data.map(obj => header.map(key => obj[key]).join(',')).join('\n');
    
        return csvContent;
    }
    
    function saveJSONToCSV(jsonData, defaultFileName) {
        // Convert JSON data to CSV
        const csvContent = JSONtoCSV(jsonData);
    
        // Create a Blob with the CSV content and download it
        const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
        const link = document.createElement('a');
        link.href = URL.createObjectURL(blob);
        link.download = defaultFileName;
        link.style.display = 'none';
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    }
    
    function reset() {
        // Reset the UI and clear stored data
    
        // Clear stored data
        protein_parameters = null;
        depolymerase_predictions = null;
        updateButtonStates();

        document.getElementById('calculate-parameters').disabled = true;
        document.getElementById('predict-depolymerases').disabled = true;

        // Reset file input elements
        document.getElementById('fasta-upload').value = '';
        document.getElementById('training-upload').value = '';
    
        // Clear the table
        $('#table-container').empty();
        
        // Clear console
        console.clear();

    }

    function updateButtonStates() {
        document.getElementById('predict-depolymerases').disabled = !protein_parameters;
        document.getElementById('save-parameters').disabled = !protein_parameters;
        document.getElementById('save-predictions').disabled = !depolymerase_predictions;
        document.getElementById('display-predictions').disabled = !depolymerase_predictions;
    }


    function setLoadingStatus(text, isComplete = false) {
        const loadingStatus = document.getElementById('loading-status');
        if (isComplete) {
            loadingStatus.innerHTML = '<span style="color: blue;">' + text + '</span>';
        } else {
            loadingStatus.innerHTML = '<div class="d-flex align-items-center"><div class="spinner-border text-primary" role="status"></div><span class="ml-2">' + text + '</span></div>';
        }
    }
    
    function disableAllButtons(disable) {
        const buttons = document.querySelectorAll('button');
        buttons.forEach((button) => {
            button.disabled = disable;
        });
    }
    
    function updateCalculateParametersButtonState() {
        const fastaFileInput = document.getElementById('fasta-upload');
        const hasFile = fastaFileInput.files.length > 0;
        document.getElementById('calculate-parameters').disabled = !hasFile;
    }

    function toggleButtonLoading(button, isLoading) {
        if (!button) {
            console.error('toggleButtonLoading: Invalid button element');
            return;
        }
        
        let spinner = button.querySelector('.spinner-border');
        
        if (isLoading) {
            if (!spinner) {
                spinner = document.createElement('span');
                spinner.className = 'spinner-border spinner-border-sm';
                spinner.setAttribute('role', 'status');
                spinner.setAttribute('aria-hidden', 'true');
                button.insertBefore(spinner, button.firstChild);
                button.dataset.originalText = button.textContent.trim();
            }
            button.disabled = true;
            button.textContent = ' Working...';
        } else {
            if (spinner) {
                button.removeChild(spinner);
            }
            button.disabled = false;
            button.textContent = button.dataset.originalText;
        }
    }

    function displayPredictionsTable() {
        const tableContainer = document.getElementById('table-container');
        const tableHTML = jsonToHTMLTable(depolymerase_predictions, 'table table-hover table-striped');
        tableContainer.innerHTML = tableHTML;
        // Initialise the Bootstrap Table plugin
        $('table').bootstrapTable();
    }
    
    function jsonToHTMLTable(jsonArray, tableClasses = '') {
        // Create container div
        const container = document.createElement('div');
    
        // Add the header (h2) element
        const header = document.createElement('h2');
        header.innerText = 'Depolymerase predictions';
        container.appendChild(header);
    
        // Create the table element
        const table = document.createElement('table');
        table.classList = tableClasses;
    
        // Add data attributes for pagination
        table.setAttribute('data-pagination', 'true');
        table.setAttribute('data-page-size', '10');
        table.setAttribute('data-side-pagination', 'client');
        table.setAttribute('data-resizable', 'true');
    
        // Table header
        const thead = document.createElement('thead');
        const headerRow = document.createElement('tr');
        if (jsonArray.length > 0) {
            const keys = Object.keys(jsonArray[0]);
            keys.forEach((key) => {
                const th = document.createElement('th');
                th.setAttribute('data-field', key);
                th.setAttribute('data-sortable', 'true');
                th.innerText = key;
                headerRow.appendChild(th);
            });
        }
        thead.appendChild(headerRow);
        table.appendChild(thead);
    
        // Table body
        const tbody = document.createElement('tbody');
        jsonArray.forEach((row) => {
            const tr = document.createElement('tr');
            for (const key in row) {
                const td = document.createElement('td');
                td.innerText = row[key];
                tr.appendChild(td);
            }
            tbody.appendChild(tr);
        });
        table.appendChild(tbody);
    
        // Append table to container div
        container.appendChild(table);
    
        return container.outerHTML;
    }


}


main();