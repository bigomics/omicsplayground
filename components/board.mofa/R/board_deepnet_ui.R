##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

DeepNetInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    ## data set parameters
    shiny::selectInput(ns("select_pheno"), "Select phenotype:", choices = NULL, multiple = FALSE),
    shiny::selectInput(ns("select_datatypes"), "Select datatypes:", choices = NULL, multiple = TRUE),
    shiny::selectInput(ns("show_datatypes"), "Show datatypes:", choices = NULL, multiple = TRUE),
    shiny::selectInput(ns("show_conditions"), "Show conditions:", choices = NULL, multiple = TRUE),
    br(),
    p("Network learning:"),
    bslib::layout_column_wrap(
      width = 1 / 2,
      gap = "5px",
      shiny::actionButton(ns("step"), "Step1", size = "xs"),
      shiny::actionButton(ns("step20"), "Step20", size = "xs"),
      shiny::actionButton(ns("step100"), "Step100", size = "xs"),
      shiny::actionButton(ns("reset"), "Reset", size = "xs")
    ),
    bslib::accordion(
      id = ns("data_type_accordion"),
      open = FALSE,
      bslib::accordion_panel(
        "Network options",
        icon = icon("cog", lib = "glyphicon"),
        shiny::tagList(
          shiny::radioButtons(ns("model"), "Model:", c("AE", "SAE", "MLP"), selected = "SAE", inline = TRUE),
          shiny::selectInput(ns("layers"), "Layers:",
            choices = c("mini", "medium", "deep"),
            selected = "mini"
          ),
          shiny::sliderInput(ns("latent_dim"), "Latent dimension:", 4, 80, 16, 8),
          shiny::checkboxInput(ns("augment"), "augment data (10x)", TRUE),
          shiny::checkboxInput(ns("addgsets"), "add genesets", FALSE),
          # shiny::checkboxInput(ns("scaleinput"), "scale input", TRUE),
          # shiny::checkboxInput(ns("sdweight"), "gradient SD weight", TRUE),
          # shiny::checkboxInput(ns("useBN"), "use BatchNorm", TRUE),
          # shiny::checkboxInput(ns("dropout"), "use dropout", FALSE),
          shiny::checkboxInput(ns("addnoise"), "add noise", TRUE),
          shiny::checkboxInput(ns("useGLU"), "use GLU", FALSE)
          # shiny::checkboxInput(ns("multitarget"), "multi target", FALSE)
          # shiny::selectInput(ns("optim"), "Optimizer",
          #  choices = c("adam","adamw","sgd","lbfgs"), selected="adam"),
          # shiny::selectInput(ns("actfun"), "Activation function",
          #  choices = c("relu","gelu","leaky","silu"), selected="leaky"),
          # shiny::sliderInput(ns("l1regularization"),"L1 regularization (log10):",-4,4,0),
          # shiny::sliderInput(ns("l2regularization"),"L2 regularization (log10):",-4,4,0)
        )
      )
    )
  )
}

DeepNetUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 700 ## full height of page
  rowH1 <- 250 ## row 1 height
  rowH2 <- 440 ## row 2 height

  shiny::div(
    boardHeader(title = "Multi-Omics Supervised Auto-Encoder", info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),

      ## ----------------------------------------------------------------
      shiny::tabPanel(
        "Model training",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>Multi-Omics integration by deep learning</b>. Here we integrate multi-omics data using a multi-view supervised auto-encoder. This architecture is a combination of (multiple) auto-encoders and a multilayer perceptron (MLP) classifier that uses a merge of the bottleneck layers for prediction.")),
          bslib::layout_columns(
            col_widths = bslib::breakpoints(
              xxxl = c(6, 6),
              xl = c(12, 12),
              sm = c(12, 12)
            ),
            height = "calc(100vh - 180px)",
            bslib::layout_columns(
              col_widths = bslib::breakpoints(
                xxxl = c(12, 6, 6),
                xl = c(6, 3, 3),
                sm = c(12, 12, 12)
              ),
              plot_deepnet_diagram_ui(
                ns("deepnet_diagram"),
                title = "Network model architecture",
                info.text = "Supervised auto-encoder (SAE) on multi-omics data performed using R torch-based deep learning pipeline. Includes feature selection, SAE models, training/validation splitting, training with distinct optimizers, prediction and later feature extraction. Depicted are the (multiple) layers of the model, including input, middle, output, and the bottleneck.",
                info.methods = "The multi-omics, log2-transformed and normalized data matrix is imputed if missing values are detected. Datatypes can be optionally selected using the {Select datatypes} option. Optionally, each data type can be 10x augmented by checking the {augment data} box. Augmentation is conducted by (i) computing average feature standard deviation; (ii) expanding the matrix (by samples) 10 times the original size; (iii) adding to the expanded matrix random noise corresponding to the product between the original average feature standard deviation and random values drawn from a Gaussian distribution. To expand the breadth of features used in the SAE, genesets from gene set enrichment analysis can be optionally added by checking the {add genesets} box as additional 'datatype'. Noise can be optionally added by checking the {add noise} box. Adding noise triggers marginal increase in feature variations without altering the biological patterns in the data and its latent space. Noise is added as the product between average feature standard deviation and random values drawn from a Gaussian distribution. Multi-omics SAE is then performed using R6 class functionalities provided in the R torch R package. A GLU activation/gating mechanism can be optionally added by checking the {use GLU} box. MultiBlockMultiTargetSAE ('MT') model is used to handle matrices of distinct features and predict multiple target variables at once. The SAE neural network is used to learn representations and make feature predictions. Internally, the model is organized into modules: (i) initialization, which converts data into torch sensors and split into training and validation; (ii) training, which sets the optimizer, defines the loss function, trains the model, and stores training and validation loss values; (iii) prediction, which runs the trained model to predict target probabilities; (iv) latent representation extraction, which returns the learned latent (encoded) representations for each view and the integrated multi-omics data; (v) feature importance by gradient, which computes feature importance for each input feature by analyzing how small changes (perturbations or gradients) may affect the output; (vi) model architecture dimensions, which returns the dimensions of encoder, decoder, and predictor layers in the model.",
                info.references = list(
                  list(
                    "Paszke, A., et al. (2019). “PyTorch: An Imperative Style, High-Performance Deep Learning Library.” arXiv:1912.01703,e1005752.", "https://doi.org/10.48550/arXiv.1912.01703"
                  ),
                  list("Le, L., Patterson, A., White, M. (2018). “Supervised autoencoders: Improving generalization performance with unsupervised regularizers.”, NeurIPS Proceedings.", "https://proceedings.neurips.cc/paper_files/paper/2018/file/2a38a4a9316c49e5a833517c45d31070-Paper.pdf")
                ),
                caption = "Supervised auto-encoder model architecture for multi-omics data. Depicted are the (multiple) layers of the model, including input, middle, output, and the bottleneck. The model is implemented using R torch-based deep learning pipeline in R. Includes feature selection, SAE models, training/validation splitting, training with distinct optimizers, prediction and later feature extraction."
              ),
              plot_deepnet_clusters_ui(
                ns("deepnet_clusters"),
                title = "Network clustering",
                info.text = "Neural net clustering. Dimensional reduction of the reduced representations ('encodings') from the multi-omics supervised auto-encoder (SAE).",
                info.methods = "Dimensional reduction of the reduced representations ('encodings') produced by the SAE model for the multi-omics data. For details on how SAE is implemented, please refer to the information provided in the network model architecture plot. Prior to dimensional reduction, small datasets are augmented by expanding the matrix (by samples) 3 times the original size. Noise is then added as the product between standard deviation and random values drawn from a Gaussian distribution. Dimensional reduction is then conducted using singular value decomposition via the R function svd. A scatter plot of the 1st component (x-axis) and 2nd component (y-axis) is presented.",
                info.references = list(
                  list(
                    "Paszke, A., et al. (2019). “PyTorch: An Imperative Style, High-Performance Deep Learning Library.” arXiv:1912.01703,e1005752.", "https://doi.org/10.48550/arXiv.1912.01703"
                  )
                ),
                caption = "Principal component analysis of the encodings from the multi-omics supervised auto-encoder (SAE). Scatter plot of the 1st principal component (x-axis) and 2nd principal component (y-axis). Each dot is a feature, coloured by its datatype."
              ),
              plot_deepnet_aescatter_ui(
                ns("deepnet_aescatter"),
                title = "Signal reconstruction",
                info.text = "Scatter plot of signal reconstruction. Reconstructed data from SAE model are compared to actual, input data to assess accuracy of SAE model.",
                info.methods = "The model's SAE output is first computed on the stored, normalized omics data views. The input data are scaled, and the SAE reconstructed data are re-combined into original dimensions. A random selection of 1000 points is performed on the scaled input and reconstructed data and a scatter plot is generated. This plot would provide a visual diagnostic of the SAE's reconstruction accuracy, enabling to assess model fit and potential under- or overfitting. It informs on how well the SAE reconstructs the input data. If the SAE performed perfectly, all points should fall along the diagonal (actual ≈ reconstructed). Deviations may indicate reconstruction error which may occur in complex biological data.",
                info.references = list(
                  list(
                    "Paszke, A., et al. (2019). “PyTorch: An Imperative Style, High-Performance Deep Learning Library.” arXiv:1912.01703,e1005752.", "https://doi.org/10.48550/arXiv.1912.01703"
                  )
                ),
                caption = "Scatter plot of actual vs SAE model-reconstructed data."
              )
            ),
            bslib::layout_columns(
              col_widths = bslib::breakpoints(
                xl = c(6, 3, 3),
                sm = c(12, 12, 12, 12)
              ),
              plot_deepnet_biomarkerheatmap_ui(
                ns("deepnet_biomarkerheatmap"),
                title = "Biomarker heatmap",
                info.text = "Heatmap of top features by (importance) across omics data types. Visualize potential biomarkers, as inferred by the SAE model, across multiple-omics data.",
                info.methods = "From the SAE model, a nested list of gradients for a phenotype and each omics view is returned. The phenotype can be selected in the option {Select phenotype}. Optionally, datatypes can be set in the {Select datatypes} option to focus on a specific datatype. For each phenotype and datatype, the features are ranked by average squared gradient (importance). The rankings across phenotypes and datatypes are then combined to select a balanced set of top 20 features. These top features are then plotted in a heatmap. For heatmap, we employ row (feature) scaling and z-score calculation, and 'ward.d2' as clustering method. This approach aids visualization of the most important features (biomarkers) across omics data types using the model’s learned gradients as feature importance scores.",
                info.references = list(
                  list(
                    "Paszke, A., et al. (2019). “PyTorch: An Imperative Style, High-Performance Deep Learning Library.” arXiv:1912.01703,e1005752.", "https://doi.org/10.48550/arXiv.1912.01703"
                  )
                ),
                caption = "Heatmap of top features by (importance) across omics data types."
              ),
              plot_deepnet_confusionmatrix_ui(
                ns("deepnet_confusionmatrix"),
                title = "Confusion matrix table",
                info.text = "The confusion matrix provides a tool to assess the performance of the classification SAE model by comparing the model’s predicted labels with the actual (ground truth) labels from the dataset.",
                info.methods = "A confusion matrix for each phenotype as selected by {Select phenotype}, is computed from the SAE-model and returned using the train or test data as selected by the plot option {Show set}. For each phenotype, the best predicted class is extracted. The true labels are then compared to the predicted labels and shown in a confusion matrix.",
                info.references = list(
                  list(
                    "Paszke, A., et al. (2019). “PyTorch: An Imperative Style, High-Performance Deep Learning Library.” arXiv:1912.01703,e1005752.", "https://doi.org/10.48550/arXiv.1912.01703"
                  )
                ),
                caption = "Table of confusion matrix. SAE model-predicted labels vs actual (ground truth) labels from the dataset."
              ),
              plot_deepnet_lossplot_ui(
                ns("deepnet_lossplot"),
                title = "Loss history",
                info.text = "Scatter plot of training and validation loss curves over the course of model training.",
                info.methods = "Scatter plot of the training and validation loss curves over the course of model training. From the SAE multi-omics model, the vector storing the training loss at each iteration (e.g., cross-entropy loss for classification), is extracted. The vector storing the validation loss at each iteration is also extracted. Zeros values are replaced with a small positive value, 10% of the minimum nonzero loss, to enable logarithm transformation. Since loss values may spans order of magnitude, the training and validation loss values are log2-transformed for better visualization. Training and validation loss values are colored in black and green, respectively. ",
                info.references = list(
                  list(
                    "Paszke, A., et al. (2019). “PyTorch: An Imperative Style, High-Performance Deep Learning Library.” arXiv:1912.01703,e1005752.", "https://doi.org/10.48550/arXiv.1912.01703"
                  )
                ),
                caption = ""
              )
            )
          )
        )
      ), ## end of tabPanel


      ## ----------------------------------------------------------------
      shiny::tabPanel(
        "Biomarker heatmap",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>Biomarker heatmap.</b>.")),
          bslib::layout_columns(
            col_widths = c(12),
            height = "calc(100vh - 180px)",
            plot_deepnet_biomarkerheatmap_ui(
              ns("deepnet_bigheatmap"),
              title = "Biomarker heatmap",
              info.text = "Heatmap of top features by (importance) across omics data types. Visualize potential biomarkers, as inferred by the SAE model, across multiple-omics data.",
              info.methods = "From the SAE model, a nested list of gradients for a phenotype and each omics view is returned. The phenotype can be selected in the option {Select phenotype}. Optionally, datatypes can be set in the {Select datatypes} option to focus on a specific datatype. For each phenotype and datatype, the features are ranked by average squared gradient (importance). The rankings across phenotypes and datatypes are then combined to select a balanced set of top features. These top features are then plotted in a heatmap. For heatmap, we employ row (feature) scaling and z-score calculation, and 'ward.d2' as clustering method. This approach aids visualization of the most important features (biomarkers) across omics data types using the model’s learned gradients as feature importance scores.",
              info.references = list(
                list(
                  "Paszke, A., et al. (2019). “PyTorch: An Imperative Style, High-Performance Deep Learning Library.” arXiv:1912.01703,e1005752.", "https://doi.org/10.48550/arXiv.1912.01703"
                )
              ),
              caption = "Heatmap of top features by (importance) across omics data types."
            )
          )
        )
      ), ## end of tabPanel

      ## ----------------------------------------------------------------
      shiny::tabPanel(
        "Gradient vs. foldchange",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>Gradient vs foldchange. </b>. This board compares input gradients of the network (i.e. change of prediction with respect to inputs) with the log-foldchange. Good biomarkers should have high foldchange and large input gradient.")),
          bslib::layout_columns(
            col_widths = c(6, 6, 12),
            height = "calc(100vh - 180px)",
            plot_deepnet_gradients_ui(
              ns("deepnet_gradients"),
              title = "Network gradients",
              info.text = "Barplot showing the most important features (by gradient magnitude) for each class in each omics view.",
              info.methods = "The gradient matrices, one for each omics view, is extracted from the SAE model. Each datatype-specific gradient matrix contains features and classes/labels for a phenotype. Classes/labels can be selected using the {Show condition} option. Phenotypes can selected using the {Select phenotype} option. Gradient values are measures of the importance of a feature in given phenotype class. A barplot showing the most important features (by gradient magnitude) for each class in each omics view is shown. Gradients are sorted for the selected features. By using the plot option {positive only}, positive and/or negative gradient values are shown.",
              info.references = list(
                list(
                  "Paszke, A., et al. (2019). “PyTorch: An Imperative Style, High-Performance Deep Learning Library.” arXiv:1912.01703,e1005752.", "https://doi.org/10.48550/arXiv.1912.01703"
                )
              ),
              caption = "Barplot of datatype-specific feature gradient values for the specified condition (class) of the selected phenotype. The datatype is indicated at the top of each barplot."
            ),
            plot_deepnet_gradients_ui(
              ns("deepnet_fcvsgrad"),
              title = "Feature's gradient vs. fold-change",
              info.text = "Scatter plot of feature fold-change vs. feature SAE gradient per datatype.",
              info.methods = "For the selected phenotype, as chosen by the {Select phenotype} option, the difference between the average expression in the phenotype classes is computed for each feature. Gene sets as inferred in the gene set enrichment analysis are also included as 'feature set'. The gradient matrices, one for each omics view, is extracted from the SAE model. For each available feature per datatype, a scatter plot of log2FC vs gradient is displayed. Feature with lowest and highest FC and gradient are coloured in blue and red, respectively.",
              # [For each datatype, the average of each squared FC and gradient feature value is computed and rank of average values computed. The two rank vectors (one for log2FC and one for gradient) are then combined by sum. This results into a feature score corresponding to the rank of log2FC and gradient. Lower scores correspond to features that have lower ranks in both FC and gradient matrices; higher scores to features that have higher log2FC and gradient ranks].
              info.references = list(
                list(
                  "Paszke, A., et al. (2019). “PyTorch: An Imperative Style, High-Performance Deep Learning Library.” arXiv:1912.01703,e1005752.", "https://doi.org/10.48550/arXiv.1912.01703"
                )
              ),
              caption = "Scatter plot of feature log2 fold-change vs. feature SAE gradient per datatype.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            table_deepnet_gradients_ui(
              ns("deepnet_table"),
              title = "Table of features' fold-change and gradient.",
              info.text = "Table of multi-omics feature log fold-change and gradient. A datatype symbols is prepended to each feature (px=proteomics; mx=metabolomics; gx=transcriptomics). Additional feature annotations are reported in the symbol and title columns. For each feature, gradient values and log2FC values for each phenotype class are reported. For details on how FC is calculated, please refer to the info of the Gradient vs. foldchange scatter plot. For details on the SAE model, please refer to the info of the 'Network model architecture' scheme under the 'Model training' tab.'",
              caption = "",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ) ## end of tabPanel
    )
  )
}
