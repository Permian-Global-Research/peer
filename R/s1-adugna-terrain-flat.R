# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
#
# """
# Version: v1$0
# Date: 2021-03-12
# Description: This code is adopted from
# Vollrath, A$, Mullissa, A$, & Reiche, J$ (2020)$
# Angular-Based Radiometric Slope Correction for Sentinel-1 on Google Earth Engine$
#   Remote Sensing, 12(11), [1867]$ https://doi$org/10$3390/rs12111867
# """
#

# ---------------------------------------------------------------------------//
# Terrain Flattening
# ---------------------------------------------------------------------------//

#' Adunga S1 slope correction
#'
#' @param collection
#' @param TERRAIN_FLATTENING_MODEL
#' @param DEM
#' @param TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER
#'
#' @return
#' @export
#'
#' @examples
slope_correction <- function(collection,
                             TERRAIN_FLATTENING_MODEL= 'VOLUME',
                             DEM=ee$Image('USGS/SRTMGL1_003'),
                             TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER=0){
  # """
  #
  #   Parameters
  #   ----------
  #   collection : ee image collection
  #       DESCRIPTION$
  #   TERRAIN_FLATTENING_MODEL : string
  #       The radiometric terrain normalization model, either volume or direct
  #   DEM : ee asset
  #       The DEM to be used
  #   TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER : integer
  #       The additional buffer to account for the passive layover and shadow
  #
  #   Returns
  #   -------
  #   ee image collection
  #       An image collection where radiometric terrain normalization is
  #       implemented on each image
  #
  #   """

  ninetyRad = ee$Image$constant(90)$multiply(pi/180)

  .volumetric_model_SCF <- function(theta_iRad, alpha_rRad){
    # """
    #
    #     Parameters
    #     ----------
    #     theta_iRad : ee$Image
    #         The scene incidence angle
    #     alpha_rRad : ee$Image
    #         Slope steepness in range
    #
    #     Returns
    #     -------
    #     ee$Image
    #         Applies the volume model in the radiometric terrain normalization
    #
    #     """

    # Volume model
    nominator = (ninetyRad$subtract(theta_iRad)$add(alpha_rRad))$tan()
    denominator = (ninetyRad$subtract(theta_iRad))$tan()
    nominator$divide(denominator)
  }


  .direct_model_SCF <- function(theta_iRad, alpha_rRad, alpha_azRad){
    # """
    #
    #     Parameters
    #     ----------
    #     theta_iRad : ee$Image
    #         The scene incidence angle
    #     alpha_rRad : ee$Image
    #         Slope steepness in range
    #
    #     Returns
    #     -------
    #     ee$Image
    #         Applies the direct model in the radiometric terrain normalization
    #
    #     """
    # Surface model
    nominator = (ninetyRad$subtract(theta_iRad))$cos()
    denominator = alpha_azRad$cos()$multiply((ninetyRad$subtract(theta_iRad)$add(alpha_rRad))$cos())
    nominator$divide(denominator)
  }


  .erode <- function(image, distance){
    # """
    #
    #
    #       Parameters
    #       ----------
    #       image : ee$Image
    #           Image to apply the erode function to
    #       distance : integer
    #           The distance to apply the buffer
    #
    #       Returns
    #       -------
    #       ee$Image
    #           An image that is masked to conpensate for passive layover
    #           and shadow depending on the given distance
    #
    #       """
    # buffer function (thanks Noel)

    d = (image$Not()$unmask(1)$fastDistanceTransform(30)$sqrt()
         $multiply(ee$Image$pixelArea()$sqrt()))

    image$updateMask(d$gt(distance))
  }


  .masking <- function(alpha_rRad, theta_iRad, buffer){
    # """
    #
    #     Parameters
    #     ----------
    #     alpha_rRad : ee$Image
    #         Slope steepness in range
    #     theta_iRad : ee$Image
    #         The scene incidence angle
    #     buffer : TYPE
    #         DESCRIPTION$
    #
    #     Returns
    #     -------
    #     ee$Image
    #         An image that is masked to conpensate for passive layover
    #         and shadow depending on the given distance
    #
    #     """
    # calculate masks
    # layover, where slope > radar viewing angle
    layover = alpha_rRad$lt(theta_iRad)$rename('layover')
    # shadow
    shadow = alpha_rRad$gt(ee$Image$constant(-1)
                           $multiply(ninetyRad$subtract(theta_iRad)))$rename('shadow')
    # combine layover and shadow
    mask = layover$And(shadow)
    # add buffer to final mask
    if (buffer > 0){
      mask = .erode(mask, buffer)
    }

    mask$rename('no_data_mask')
  }



  .correct <- function(image){
    # """
    #
    #
    #     Parameters
    #     ----------
    #     image : ee$Image
    #         Image to apply the radiometric terrain normalization to
    #
    #     Returns
    #     -------
    #     ee$Image
    #         Radiometrically terrain corrected image
    #
    #     """

    bandNames = image$bandNames()

    geom = image$geometry()
    proj = image$select(1)$projection()

    elevation = DEM$resample('bilinear')$reproject(proj,NULL, 10)$clip(geom)

    # calculate the look direction
    heading = ee$Terrain$aspect(image$select('angle'))$reduceRegion(ee$Reducer$mean(), image$geometry(), 1000)


    #in case of null values for heading replace with 0
    heading = ee$Dictionary(heading)$combine(list('aspect'= 0), FALSE)$get('aspect')

    heading = ee$Algorithms$If(
      ee$Number(heading)$gt(180),
      ee$Number(heading)$subtract(360),
      ee$Number(heading)
    )

    # the numbering follows the article chapters
    # 2$1$1 Radar geometry
    theta_iRad = image$select('angle')$multiply(pi/180)
    phi_iRad = ee$Image$constant(heading)$multiply(pi/180)

    # 2$1$2 Terrain geometry
    alpha_sRad = ee$Terrain$slope(elevation)$select('slope')$multiply(pi / 180)

    aspect = ee$Terrain$aspect(elevation)$select('aspect')$clip(geom)

    aspect_minus = aspect$updateMask(aspect$gt(180))$subtract(360)

    phi_sRad = aspect$updateMask(aspect$lte(180))$
      unmask()$
      add(aspect_minus$unmask())$
      multiply(-1)$
      multiply(pi / 180)

    #elevation = DEM$reproject(proj,NULL, 10)$clip(geom)

    # 2$1$3 Model geometry
    # reduce to 3 angle
    phi_rRad = phi_iRad$subtract(phi_sRad)

    # slope steepness in range (eq$ 2)
    alpha_rRad = (alpha_sRad$tan()$multiply(phi_rRad$cos()))$atan()

    # slope steepness in azimuth (eq 3)
    alpha_azRad = (alpha_sRad$tan()$multiply(phi_rRad$sin()))$atan()

    # 2$2
    # Gamma_nought
    gamma0 = image$divide(theta_iRad$cos())

    if (TERRAIN_FLATTENING_MODEL == 'VOLUME'){
      # Volumetric Model
      scf = .volumetric_model_SCF(theta_iRad, alpha_rRad)
    } else if (TERRAIN_FLATTENING_MODEL == 'DIRECT'){
      scf = .direct_model_SCF(theta_iRad, alpha_rRad, alpha_azRad)
    }


    # apply model for Gamm0
    gamma0_flat = gamma0$multiply(scf)

    # get Layover/Shadow mask
    mask = .masking(alpha_rRad, theta_iRad, TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER)
    output = gamma0_flat$mask(mask)$rename(bandNames)$copyProperties(image)
    output = ee$Image(output)$addBands(image$select('angle'), NULL, TRUE)

    output$set('system:time_start', image$get('system:time_start'))
  }

 collection$map(.correct)
 # .correct(collection)
}
