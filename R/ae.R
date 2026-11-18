# Autoencoder in base R: tanh hidden layers, linear output, MSE loss,
# analytic backpropagation and Adam. dims = c(d, k, d) gives the shallow
# model; c(d, 8, k, 8, d) gives the deeper variant.

init_weights <- function(din, dout) {
  limit <- sqrt(6 / (din + dout))
  matrix(runif(din * dout, -limit, limit), nrow = din, ncol = dout)
}

ae_init <- function(dims) {
  purrr::map(seq_len(length(dims) - 1), \(l) {
    list(W = init_weights(dims[l], dims[l + 1]), b = rep(0, dims[l + 1]))
  })
}

ae_forward <- function(params, X) {
  n_layers <- length(params)
  acts <- vector("list", n_layers + 1)
  acts[[1]] <- X
  for (l in seq_len(n_layers)) {
    H <- acts[[l]] %*% params[[l]]$W + rep(params[[l]]$b, each = nrow(X))
    acts[[l + 1]] <- if (l < n_layers) tanh(H) else H
  }
  acts
}

# Loss = (1 / (N d)) * ||X - X_hat||_F^2 and its gradients
ae_loss_grad <- function(params, X) {
  acts <- ae_forward(params, X)
  n_layers <- length(params)
  diff <- acts[[n_layers + 1]] - X
  loss <- sum(diff^2) / length(X)
  delta <- 2 * diff / length(X)
  grads <- vector("list", n_layers)
  for (l in n_layers:1) {
    grads[[l]] <- list(W = crossprod(acts[[l]], delta), b = colSums(delta))
    if (l > 1) delta <- tcrossprod(delta, params[[l]]$W) * (1 - acts[[l]]^2)
  }
  list(loss = loss, grads = grads)
}

ae_mse <- function(params, X) {
  acts <- ae_forward(params, X)
  mean((acts[[length(acts)]] - X)^2)
}

train_ae <- function(X_train, X_val, dims, max_epochs = 60, patience = 5,
                     batch_size = 256, lr = 0.005, beta1 = 0.9, beta2 = 0.999,
                     eps = 1e-8) {
  params <- ae_init(dims)
  zeros <- purrr::map(params, \(p) list(W = p$W * 0, b = p$b * 0))
  m <- zeros
  v <- zeros
  step <- 0
  N <- nrow(X_train)
  n_batches <- ceiling(N / batch_size)
  train_hist <- val_hist <- rep(NA_real_, max_epochs)
  best <- list(loss = Inf, params = params, epoch = 0L)

  for (epoch in seq_len(max_epochs)) {
    idx <- sample.int(N)
    batch_losses <- numeric(n_batches)
    for (b in seq_len(n_batches)) {
      rows <- idx[((b - 1) * batch_size + 1):min(b * batch_size, N)]
      lg <- ae_loss_grad(params, X_train[rows, , drop = FALSE])
      batch_losses[b] <- lg$loss
      step <- step + 1
      for (l in seq_along(params)) {
        for (nm in c("W", "b")) {
          g <- lg$grads[[l]][[nm]]
          m[[l]][[nm]] <- beta1 * m[[l]][[nm]] + (1 - beta1) * g
          v[[l]][[nm]] <- beta2 * v[[l]][[nm]] + (1 - beta2) * g^2
          m_hat <- m[[l]][[nm]] / (1 - beta1^step)
          v_hat <- v[[l]][[nm]] / (1 - beta2^step)
          params[[l]][[nm]] <- params[[l]][[nm]] - lr * m_hat / (sqrt(v_hat) + eps)
        }
      }
    }
    # Training loss averages mini-batches while weights move; validation loss
    # is computed once at the end of the epoch.
    train_hist[epoch] <- mean(batch_losses)
    val_hist[epoch] <- ae_mse(params, X_val)
    if (val_hist[epoch] < best$loss) {
      best <- list(loss = val_hist[epoch], params = params, epoch = epoch)
    } else if (epoch - best$epoch >= patience) {
      break
    }
  }

  list(
    params = best$params, dims = dims, best_epoch = best$epoch,
    epochs_run = epoch, best_val_loss = best$loss,
    train_losses = train_hist[seq_len(epoch)],
    val_losses = val_hist[seq_len(epoch)]
  )
}

# Per-sample loss = (1/d) * ||x - x_hat||^2
ae_reconstruct <- function(model, X) {
  acts <- ae_forward(model$params, X)
  diff <- acts[[length(acts)]] - X
  list(residuals = abs(diff), sample_mse = rowMeans(diff^2))
}

# Central finite-difference check of every gradient block
gradient_check <- function(params, X, h = 1e-5) {
  analytic <- ae_loss_grad(params, X)$grads
  purrr::map(seq_along(params), \(l) {
    purrr::map(c("W", "b"), \(nm) {
      theta <- params[[l]][[nm]]
      num <- theta * 0
      for (i in seq_along(theta)) {
        p_plus <- params
        p_minus <- params
        p_plus[[l]][[nm]][i] <- theta[i] + h
        p_minus[[l]][[nm]][i] <- theta[i] - h
        num[i] <- (ae_mse(p_plus, X) - ae_mse(p_minus, X)) / (2 * h)
      }
      a <- analytic[[l]][[nm]]
      tibble::tibble(
        layer = l, block = nm, n_params = length(theta),
        rel_error = sqrt(sum((a - num)^2)) / (sqrt(sum(a^2)) + sqrt(sum(num^2)))
      )
    }) |>
      dplyr::bind_rows()
  }) |>
    dplyr::bind_rows()
}
