#!/usr/bin/env python

"""Variational autoencoder model
"""

import numpy as np
from tensorflow.keras import backend as K
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Dense, Lambda
from tensorflow.keras.losses import binary_crossentropy
from tensorflow.keras.optimizers import Adam


def sampling(args):
    """Sample from normal distribution.
    """
    z_mean, z_log_var, latent_dim = args
    epsilon = K.random_normal(
        shape=(K.shape(z_mean)[0], latent_dim), mean=0., stddev=1.
    )
    return z_mean + K.exp(z_log_var) * epsilon


def vae_encoder(input_dim, neuron_nums, latent_dim=2):
    """Build VAE encoder model.

    Parameters
    ----------
    input_dim : tuple
        Input space dimensions.
    neuron_nums : List[int]
        Number of neurons in each layer.
    latent_dim : int (optional)
        Latent space dimensions.

    Returns
    -------
    encoder : keras.Model
        Encoder model.
    z_mean : keras.Model.Dense
    z_log_var : keras.Model.Dense
    encoder_input : keras.Model.Input
        Encoder input layer
    """
    encoder_input = layer = Input(shape=input_dim)

    for neuron_num in neuron_nums:
        layer = Dense(neuron_num, activation='relu')(layer)

    z_mean = Dense(latent_dim)(layer)
    z_log_var = Dense(latent_dim)(layer)

    z = Lambda(sampling)([z_mean, z_log_var, latent_dim])
    encoder = Model(encoder_input, [z_mean, z_log_var, z])

    return encoder, z_mean, z_log_var, encoder_input


def vae_decoder(input_dim, neuron_nums, latent_dim=2):
    """Build VAE decoder model.

    Parameters
    ----------
    input_dim : tuple
        Input space dimensions.
    neuron_nums : List[int]
        Number of neurons in each layer.
    latent_dim : int (optional)
        Latent space dimensions.

    Returns
    --------
    decoder : keras.Model
        Decoder model.
    """
    latent_inputs = layer = Input(shape=(latent_dim,))

    for neuron_num in neuron_nums:
        layer = Dense(neuron_num, activation='relu')(layer)

    layer = Dense(np.prod(input_dim), activation='sigmoid')(layer)

    decoder = Model(latent_inputs, layer, name='decoder')
    return decoder


def build_vae(
    input_dim, encoder_neuron_nums,
    decoder_neuron_nums=None, latent_dim=2
):
    """Build VAE model.

    Parameters
    ----------
    input_dim : tuple
        Input space dimensions.
    encoder_neuron_nums : List[int]
        Number of neurons in encoder module.
    decoder_neuron_nums : List[int] (optional)
        Number of neurons in decoder module.
        By default, it's the reverse of encoder_neuron_nums.
    latent_dim : int (optional)
        Latent space dimensions.

    Returns
    -------
    encoder : keras.Model
        constructed encoder model
    decoder : keras.Model
        constructed decoder model
    vae : keras.Model
        constructed VAE model
    """
    encoder, z_mean, z_log_var, encoder_input = vae_encoder(
        input_dim, encoder_neuron_nums, latent_dim)

    # if decoder_neuron_nums is not given,
    # set to the reverse of encoder_neuron_nums
    if decoder_neuron_nums is None:
        decoder_neuron_nums = encoder_neuron_nums[::-1]

    decoder = vae_decoder(input_dim, decoder_neuron_nums, latent_dim)

    z_decoded = decoder(encoder(encoder_input)[2])

    def vae_loss(encoder_input, z_decoded):
        # reconstruction loss
        xent_loss = binary_crossentropy(encoder_input, z_decoded)
        xent_loss *= np.prod(input_dim)

        kl_loss = -0.5 * K.sum(
            1 + z_log_var - K.square(z_mean) - K.exp(z_log_var), axis=-1
        )

        return K.mean(xent_loss + kl_loss)

    vae = Model(encoder_input, z_decoded)
    vae.add_loss(vae_loss(encoder_input, z_decoded))
    vae.compile(optimizer=Adam())

    return encoder, decoder, vae
